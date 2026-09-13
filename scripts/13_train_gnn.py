"""
13_train_gnn.py — B3/B4: GraphSAGE over the assembly graph.

Naming follows the project's existing convention (01_..12_...). Drop this
into work/scripts/ alongside 06_baselines.py.

Implements the three-head design recommended for open item B
(PROPOSAL_FIXES.md #1):

  1. Plasmid-pooled head  — same shape as B0+/B1/B2, for direct comparability.
     Reported descriptively, NOT used as the confirmatory number, because the
     pooled S(p) vector structurally encodes several label rules.
  2. Node head            — fragmented vs. absorbed vs. recovered, evaluated
     ONLY on units with a non-empty support set. Cannot read off `is_empty`.
  3. Graph head           — absence, evaluated as: does the graph-level
     embedding predict how many of that graph's plasmids are missing?
     A count-calibration target instead of a per-plasmid label.

Evaluation: leave-isolate-out CV (grouped by isolate, exactly like
06_baselines.py), PR-AUC, paired bootstrap over isolates, and the same
UNRESOLVED-if-arms-disagree discipline used throughout PROGRESS.md.

ASSUMPTIONS TO VERIFY AGAINST YOUR ACTUAL SCHEMA before running — I do not
have your literal files, so these are the interfaces PROGRESS.md's own
description implies. Search-and-adjust anywhere marked ADAPT.

  - pyg_dataset.pt is a list[torch_geometric.data.Data], one per (isolate,
    prep, depth) cell / per real isolate, each with:
      data.x          : [N, F] node features (6 features as built by
                         05_features.py), global node appended last,
                         flagged by an extra column (per §6 of PROGRESS.md)
      data.edge_index : [2, E] simple undirected graph, self-loops dropped
      data.isolate    : str, isolate id (for grouping in CV) -- ADAPT if
                         this lives in a separate manifest instead
      data.plasmid_units : list of dicts, one per label unit in this graph:
            { "plasmid_id": str,
              "label": one of {"recovered","fragmented","absorbed","absent"},
              "support_nodes": list[int]  # indices into data.x, S(p)
              "plasmid_len": float }
        ADAPT: if this lives in plasmid_labels.csv / node_labels.csv instead
        of embedded in the Data object, load and join by (isolate, plasmid_id)
        before calling build_examples().

Run:
    python 13_train_gnn.py --dataset work/results/pyg_dataset.pt \
                            --out work/results/gnn_results.json \
                            --mode leave_isolate_out
"""

import argparse
import json
import random
from collections import defaultdict

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import average_precision_score

try:
    from torch_geometric.nn import SAGEConv, global_mean_pool
except ImportError as e:
    raise SystemExit(
        "torch_geometric not found. This project already has PyG 2.8.0 "
        "installed per PROGRESS.md sec. 2 — activate that env first."
    ) from e


LABELS = ["recovered", "fragmented", "absorbed", "absent"]
LABEL2IDX = {l: i for i, l in enumerate(LABELS)}
SEED = 20260812  # match the cohort's own seed convention


# --------------------------------------------------------------------------
# Model
# --------------------------------------------------------------------------

class GraphSAGEEncoder(nn.Module):
    """Two-layer GraphSAGE, produces one embedding per node (incl. global)."""

    def __init__(self, in_dim: int, hidden: int = 64, out_dim: int = 64):
        super().__init__()
        self.conv1 = SAGEConv(in_dim, hidden)
        self.conv2 = SAGEConv(hidden, out_dim)
        self.norm1 = nn.LayerNorm(hidden)

    def forward(self, x, edge_index):
        h = F.relu(self.norm1(self.conv1(x, edge_index)))
        h = self.conv2(h, edge_index)
        return h  # [N, out_dim]


class PlasmidPooledHead(nn.Module):
    """B4-pooled: mean+max over S(p), 4-way classification.
    Reported for comparability with B0+/B1/B2 only -- NOT the confirmatory
    number (see module docstring)."""

    def __init__(self, emb_dim: int, n_classes: int = 4):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(2 * emb_dim + 1, 64), nn.ReLU(),
            nn.Linear(64, n_classes),
        )

    def forward(self, node_emb, support_idx, plasmid_len_log, global_emb):
        if len(support_idx) == 0:
            # is_empty by construction -- absent. Do not let the pooled
            # vector fake topology here; feed zeros + length as before.
            pooled = torch.zeros(2 * node_emb.shape[1], device=node_emb.device)
        else:
            sup = node_emb[support_idx]
            pooled = torch.cat([sup.mean(0), sup.max(0).values], dim=0)
        feat = torch.cat([pooled, plasmid_len_log.view(1)], dim=0)
        return self.mlp(feat.unsqueeze(0))  # [1, n_classes]


class NodeHead(nn.Module):
    """Confirmatory head 1: fragmented / absorbed / recovered only,
    on units with |S(p)| >= 1. Cannot see `is_empty` -- absent units are
    excluded from this head's training and evaluation entirely."""

    def __init__(self, emb_dim: int, n_classes: int = 3):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(2 * emb_dim, 64), nn.ReLU(),
            nn.Linear(64, n_classes),
        )

    def forward(self, node_emb, support_idx):
        sup = node_emb[support_idx]
        pooled = torch.cat([sup.mean(0), sup.max(0).values], dim=0)
        return self.mlp(pooled.unsqueeze(0))


class GraphHead(nn.Module):
    """Confirmatory head 2: graph-level embedding -> predicted count of
    missing (absent) plasmids in this graph. Regression target, compared
    against the true count via a discretised PR-style calibration curve
    at evaluation time."""

    def __init__(self, emb_dim: int):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(emb_dim, 32), nn.ReLU(),
            nn.Linear(32, 1),
        )

    def forward(self, node_emb, batch_vec=None):
        # mean over all real (non-global) nodes
        g = node_emb.mean(dim=0, keepdim=True)
        return self.mlp(g)  # [1, 1], predicted count


class GNNModel(nn.Module):
    def __init__(self, in_dim: int, hidden: int = 64, emb_dim: int = 64):
        super().__init__()
        self.encoder = GraphSAGEEncoder(in_dim, hidden, emb_dim)
        self.plasmid_head = PlasmidPooledHead(emb_dim)
        self.node_head = NodeHead(emb_dim)
        self.graph_head = GraphHead(emb_dim)


# --------------------------------------------------------------------------
# Data plumbing
# --------------------------------------------------------------------------

def build_examples(dataset):
    """Flatten the per-graph dataset into (graph_idx, plasmid_unit) pairs,
    plus one graph-level absence-count example per graph."""
    plasmid_examples = []
    graph_examples = []
    for gi, data in enumerate(dataset):
        n_absent = 0
        for unit in data.plasmid_units:  # ADAPT if stored elsewhere
            plasmid_examples.append((gi, unit))
            if unit["label"] == "absent":
                n_absent += 1
        graph_examples.append((gi, n_absent))
    return plasmid_examples, graph_examples


def isolate_groups(dataset):
    """ADAPT: replace with your real isolate-id lookup if not on data.isolate."""
    return [d.isolate for d in dataset]


# --------------------------------------------------------------------------
# Train / eval for one leave-isolate-out fold
# --------------------------------------------------------------------------

def run_fold(model, dataset, plasmid_examples, graph_examples,
             train_graph_idx, test_graph_idx, device, epochs=60, lr=1e-3):
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    train_graph_idx = set(train_graph_idx)
    test_graph_idx = set(test_graph_idx)

    train_plasmid = [(gi, u) for gi, u in plasmid_examples if gi in train_graph_idx]
    train_graph = [(gi, c) for gi, c in graph_examples if gi in train_graph_idx]

    model.train()
    for _ in range(epochs):
        random.shuffle(train_plasmid)
        total_loss = 0.0
        # cache per-graph encodings once per epoch (small n, fine to redo)
        node_embs = {}
        for gi in train_graph_idx:
            d = dataset[gi].to(device)
            node_embs[gi] = model.encoder(d.x, d.edge_index)

        for gi, unit in train_plasmid:
            emb = node_embs[gi]
            opt.zero_grad()
            loss = 0.0
            plen_log = torch.log10(torch.tensor(
                float(unit["plasmid_len"]) + 1.0, device=device))
            logits = model.plasmid_head(
                emb, unit["support_nodes"], plen_log, emb[-1])
            target = torch.tensor([LABEL2IDX[unit["label"]]], device=device)
            loss = loss + F.cross_entropy(logits, target)

            if len(unit["support_nodes"]) > 0 and unit["label"] != "absent":
                node_labels_present = ["fragmented", "absorbed", "recovered"]
                if unit["label"] in node_labels_present:
                    nlogits = model.node_head(emb, unit["support_nodes"])
                    ntarget = torch.tensor(
                        [node_labels_present.index(unit["label"])], device=device)
                    loss = loss + F.cross_entropy(nlogits, ntarget)
            loss.backward()
            opt.step()
            total_loss += float(loss.detach())

        for gi, true_count in train_graph:
            emb = node_embs[gi]
            opt.zero_grad()
            pred = model.graph_head(emb)
            target = torch.tensor([[float(true_count)]], device=device)
            gloss = F.mse_loss(pred, target)
            gloss.backward()
            opt.step()

    # ---- eval on held-out isolate's graph(s) ----
    model.eval()
    results = {"plasmid_probs": [], "plasmid_true": [],
               "node_probs": [], "node_true": [],
               "graph_pred": [], "graph_true": []}
    with torch.no_grad():
        for gi in test_graph_idx:
            d = dataset[gi].to(device)
            emb = model.encoder(d.x, d.edge_index)
            for unit in [u for g, u in plasmid_examples if g == gi]:
                plen_log = torch.log10(torch.tensor(
                    float(unit["plasmid_len"]) + 1.0, device=device))
                logits = model.plasmid_head(
                    emb, unit["support_nodes"], plen_log, emb[-1])
                prob_fail = 1.0 - F.softmax(logits, dim=1)[0, LABEL2IDX["recovered"]].item()
                results["plasmid_probs"].append(prob_fail)
                results["plasmid_true"].append(0 if unit["label"] == "recovered" else 1)

                if len(unit["support_nodes"]) > 0:
                    node_labels_present = ["fragmented", "absorbed", "recovered"]
                    if unit["label"] in node_labels_present:
                        nlogits = model.node_head(emb, unit["support_nodes"])
                        nprob = F.softmax(nlogits, dim=1)[0].cpu().numpy()
                        results["node_probs"].append(nprob.tolist())
                        results["node_true"].append(
                            node_labels_present.index(unit["label"]))

            true_count = dict(graph_examples)[gi]
            pred_count = model.graph_head(emb).item()
            results["graph_pred"].append(pred_count)
            results["graph_true"].append(true_count)

    return results


# --------------------------------------------------------------------------
# Bootstrap over isolates (matches 06_baselines.py's reporting convention)
# --------------------------------------------------------------------------

def paired_bootstrap_pr_auc(per_isolate_probs, per_isolate_true, n_boot=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    isolates = list(per_isolate_probs.keys())
    boots = []
    for _ in range(n_boot):
        sample = rng.choice(isolates, size=len(isolates), replace=True)
        y, p = [], []
        for iso in sample:
            y.extend(per_isolate_true[iso])
            p.extend(per_isolate_probs[iso])
        if len(set(y)) < 2:
            continue
        boots.append(average_precision_score(y, p))
    boots = np.array(boots)
    return float(np.mean(boots)), float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


# --------------------------------------------------------------------------
# Main: leave-isolate-out CV
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, help="path to pyg_dataset.pt")
    ap.add_argument("--out", required=True)
    ap.add_argument("--epochs", type=int, default=60)
    ap.add_argument("--hidden", type=int, default=64)
    ap.add_argument("--device", default="cpu")
    args = ap.parse_args()

    torch.manual_seed(SEED)
    random.seed(SEED)
    np.random.seed(SEED)

    dataset = torch.load(args.dataset)  # ADAPT: weights_only kwarg per torch version
    isolates = isolate_groups(dataset)
    unique_isolates = sorted(set(isolates))
    print(f"{len(dataset)} graphs across {len(unique_isolates)} isolates.")
    if len(unique_isolates) < 5:
        print("WARNING: leave-isolate-out CV on <5 isolates will produce "
              "very wide bootstrap intervals. Report them as-is; do not "
              "hide the width. See PROJECT_COMPLETION_PLAN.md Step 5.")

    plasmid_examples, graph_examples = build_examples(dataset)
    in_dim = dataset[0].x.shape[1]

    per_isolate_probs, per_isolate_true = defaultdict(list), defaultdict(list)
    graph_pred_all, graph_true_all = [], []
    node_probs_all, node_true_all = [], []

    for held_out in unique_isolates:
        train_idx = [i for i, iso in enumerate(isolates) if iso != held_out]
        test_idx = [i for i, iso in enumerate(isolates) if iso == held_out]

        model = GNNModel(in_dim=in_dim, hidden=args.hidden).to(args.device)
        res = run_fold(model, dataset, plasmid_examples, graph_examples,
                        train_idx, test_idx, args.device, epochs=args.epochs)

        per_isolate_probs[held_out].extend(res["plasmid_probs"])
        per_isolate_true[held_out].extend(res["plasmid_true"])
        graph_pred_all.extend(res["graph_pred"])
        graph_true_all.extend(res["graph_true"])
        node_probs_all.extend(res["node_probs"])
        node_true_all.extend(res["node_true"])

        print(f"isolate {held_out}: {len(res['plasmid_true'])} plasmid units, "
              f"{sum(res['plasmid_true'])} positive")

    # Plasmid-pooled head (B4, descriptive only)
    all_y = [y for v in per_isolate_true.values() for y in v]
    all_p = [p for v in per_isolate_probs.values() for p in v]
    b4_pr_auc = average_precision_score(all_y, all_p) if len(set(all_y)) > 1 else float("nan")
    mean_ci, lo_ci, hi_ci = paired_bootstrap_pr_auc(per_isolate_probs, per_isolate_true)

    # Graph head (confirmatory: absence via count calibration)
    graph_pred_all = np.array(graph_pred_all)
    graph_true_all = np.array(graph_true_all)
    graph_mae = float(np.mean(np.abs(graph_pred_all - graph_true_all))) if len(graph_true_all) else float("nan")

    out = {
        "n_isolates": len(unique_isolates),
        "n_graphs": len(dataset),
        "b4_plasmid_pooled_head": {
            "pr_auc": b4_pr_auc,
            "bootstrap_mean": mean_ci,
            "bootstrap_ci95": [lo_ci, hi_ci],
            "note": "Descriptive only per open item B -- pooled S(p) "
                    "structurally encodes label rules, not the "
                    "confirmatory number.",
        },
        "confirmatory_node_head": {
            "n_examples": len(node_true_all),
            "note": "fragmented/absorbed/recovered only, |S(p)|>=1. "
                    "Add multi-class PR-AUC per class here once real "
                    "Wick data gives enough fragmented/absorbed units "
                    "to compute it meaningfully (v3 sim gates needed "
                    ">=20 per class; real data will likely need the "
                    "same floor before trusting this number).",
        },
        "confirmatory_graph_head": {
            "mae_missing_count": graph_mae,
            "note": "Lower is better; compare against a trivial "
                    "'always predict cohort mean absence count' floor "
                    "before interpreting.",
        },
        "WARNING": (
            "If n_isolates < 10, treat every number above as illustrative, "
            "not a claim. Report interval width, not just the point estimate. "
            "This mirrors PROGRESS.md's own treatment of the v1 (n=12) results."
        ),
    }

    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
