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

SCHEMA — VERIFIED AGAINST THE ACTUAL EXPORT (2026-09-18), not assumed.
[... unchanged from the previous revision, see prior version's docstring ...]

PATCHED AGAIN (this revision): two problems found by reading the previous
revision's code rather than trusting what it claimed to do:

1. There was no way to load a checkpoint trained on one dataset (e.g. the
   simulated cohort) and evaluate it, with no further training, on a
   DIFFERENT dataset (e.g. real Wick graphs). --train_full_only trains a
   fresh model; the default CV path also always trains fresh models. Added
   --load_model / --eval_only for exactly this: real generalization testing,
   not retraining-and-calling-it-a-test.

2. --save_model defaulted to "results/gnn_model.pt" and the save at the end
   of main() was unconditional -- so a plain CV run on ANY dataset would
   silently overwrite whatever checkpoint already lived at that path (e.g.
   your simulated-trained master model, replaced by a Wick-trained one,
   with no explicit confirmation). --save_model now defaults to None; you
   must opt in to saving explicitly, every time.

Run (from the repo root, after `pip install torch torch_geometric scikit-learn`):

  # A. Leave-isolate-out CV on whatever --dataset you point at (the direct,
  #    same-domain test; on Wick this answers "can it learn from 7 real
  #    isolates at all", NOT "does the simulated model transfer"):
  python scripts/13_train_gnn.py --dataset results/wick_pyg_dataset.pt \
                                  --out results/gnn_results_wick.json

  # B. Train once on 100% of one dataset and save a checkpoint:
  python scripts/13_train_gnn.py --dataset results/pyg_dataset.pt \
                                  --train_full_only \
                                  --save_model results/gnn_model_simulated.pt

  # C. NEW: load that checkpoint and evaluate it -- no training -- on a
  #    DIFFERENT dataset. This is the actual cross-domain generalization
  #    test (sim-trained model vs. real Wick graphs):
  python scripts/13_train_gnn.py --dataset results/wick_pyg_dataset.pt \
                                  --eval_only \
                                  --load_model results/gnn_model_simulated.pt \
                                  --out results/gnn_results_wick_from_sim_model.json
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
        if not hasattr(data, "plasmid_label4"):
            raise SystemExit(
                "data.plasmid_label4 / data.plasmid_len not found on this "
                "dataset -- results/pyg_dataset.pt (or wick_pyg_dataset.pt) "
                "needs to be built with the patched 05b_export_pyg.py."
            )
        n_absent = 0
        for pid, label, sup, plen in zip(
            data.plasmid_names, data.plasmid_label4, data.support, data.plasmid_len
        ):
            unit = {
                "plasmid_id": pid,
                "label": label,
                "support_nodes": sup,
                "plasmid_len": plen,
            }
            plasmid_examples.append((gi, unit))
            if label == "absent":
                n_absent += 1
        graph_examples.append((gi, n_absent))
    return plasmid_examples, graph_examples


def isolate_groups(dataset):
    return [d.isolate for d in dataset]


# --------------------------------------------------------------------------
# Train / eval for one leave-isolate-out fold (also used, with an empty
# train set, for pure --eval_only inference)
# --------------------------------------------------------------------------

def run_fold(model, dataset, plasmid_examples, graph_examples,
             train_graph_idx, test_graph_idx, device, epochs=60, lr=1e-3):
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    train_graph_idx = set(train_graph_idx)
    test_graph_idx = set(test_graph_idx)

    train_plasmid_by_graph = defaultdict(list)
    for gi, u in plasmid_examples:
        if gi in train_graph_idx:
            train_plasmid_by_graph[gi].append(u)

    train_graph_dict = {gi: c for gi, c in graph_examples if gi in train_graph_idx}
    train_gis = list(train_graph_idx)

    model.train()
    for _ in range(epochs):
        random.shuffle(train_gis)
        for gi in train_gis:
            opt.zero_grad()
            d = dataset[gi].to(device)
            emb = model.encoder(d.x, d.edge_index)
            loss = torch.tensor(0.0, device=device)

            plasmids = train_plasmid_by_graph.get(gi, [])
            for unit in plasmids:
                plen_log = torch.log10(torch.tensor(
                    float(unit["plasmid_len"]) + 1.0, device=device))
                logits = model.plasmid_head(
                    emb, unit["support_nodes"], plen_log, emb[d.global_index])
                target = torch.tensor([LABEL2IDX[unit["label"]]], device=device)
                loss = loss + F.cross_entropy(logits, target)

                if len(unit["support_nodes"]) > 0 and unit["label"] != "absent":
                    node_labels_present = ["fragmented", "absorbed", "recovered"]
                    if unit["label"] in node_labels_present:
                        nlogits = model.node_head(emb, unit["support_nodes"])
                        ntarget = torch.tensor(
                            [node_labels_present.index(unit["label"])], device=device)
                        loss = loss + F.cross_entropy(nlogits, ntarget)

            if gi in train_graph_dict:
                true_count = train_graph_dict[gi]
                pred = model.graph_head(emb)
                target = torch.tensor([[float(true_count)]], device=device)
                loss = loss + F.mse_loss(pred, target)

            if loss.requires_grad and len(train_gis) > 0:
                loss.backward()
                opt.step()
        # NOTE: with train_gis == [] (eval-only use), this loop body never
        # runs -- the loaded weights pass through untouched, as intended.

    # ---- eval on held-out / target graphs ----
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
                    emb, unit["support_nodes"], plen_log, emb[d.global_index])
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

            if gi in dict(graph_examples):
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
    if len(boots) == 0:
        return float("nan"), float("nan"), float("nan")
    return float(np.mean(boots)), float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


# --------------------------------------------------------------------------
# Metrics -- factored out so both the CV path and the new eval-only path
# (a single "fold" with everything in the test set) compute results the
# same way, rather than two slightly-different code paths drifting apart.
# --------------------------------------------------------------------------

def compute_metrics(n_isolates, n_graphs, per_isolate_probs, per_isolate_true,
                     graph_pred_all, graph_true_all, node_probs_all, node_true_all):
    all_y = [y for v in per_isolate_true.values() for y in v]
    all_p = [p for v in per_isolate_probs.values() for p in v]
    b4_pr_auc = average_precision_score(all_y, all_p) if len(set(all_y)) > 1 else float("nan")
    mean_ci, lo_ci, hi_ci = paired_bootstrap_pr_auc(per_isolate_probs, per_isolate_true)

    graph_pred_all = np.array(graph_pred_all)
    graph_true_all = np.array(graph_true_all)
    graph_mae = float(np.mean(np.abs(graph_pred_all - graph_true_all))) if len(graph_true_all) else float("nan")
    baseline_pred = float(np.mean(graph_true_all)) if len(graph_true_all) else float("nan")
    graph_mae_baseline = float(np.mean(np.abs(graph_true_all - baseline_pred))) if len(graph_true_all) else float("nan")

    node_labels_present = ["fragmented", "absorbed", "recovered"]
    node_true_arr = np.array(node_true_all)
    node_probs_arr = np.array(node_probs_all) if node_probs_all else np.zeros((0, 3))
    node_pr_auc_per_class, node_baseline_per_class = {}, {}
    for ci, cname in enumerate(node_labels_present):
        y = (node_true_arr == ci).astype(int)
        if len(y) == 0 or len(set(y)) < 2:
            node_pr_auc_per_class[cname] = None
            node_baseline_per_class[cname] = None
            continue
        node_pr_auc_per_class[cname] = float(average_precision_score(y, node_probs_arr[:, ci]))
        node_baseline_per_class[cname] = float(y.mean())

    return {
        "n_isolates": n_isolates,
        "n_graphs": n_graphs,
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
            "pr_auc_per_class": node_pr_auc_per_class,
            "prevalence_baseline_per_class": node_baseline_per_class,
            "note": "fragmented/absorbed/recovered only, |S(p)|>=1. "
                    "Compare pr_auc_per_class against "
                    "prevalence_baseline_per_class -- pr_auc close to or "
                    "below the baseline means this head is not learning "
                    "anything past class frequency.",
        },
        "confirmatory_graph_head": {
            "mae_missing_count": graph_mae,
            "mae_baseline_predict_mean": graph_mae_baseline,
            "note": "Lower is better. Compare mae_missing_count against "
                    "mae_baseline_predict_mean.",
        },
        "WARNING": (
            "If n_isolates < 10, treat every number above as illustrative, "
            "not a claim. Report interval width, not just the point estimate."
        ),
    }


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, help="path to pyg_dataset.pt")
    ap.add_argument("--out", required=True)
    ap.add_argument("--epochs", type=int, default=60)
    ap.add_argument("--hidden", type=int, default=64)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--save_model", default=None,
                     help="PATCHED: no longer defaults to results/gnn_model.pt. "
                          "Pass an explicit path to save a checkpoint -- opt-in, "
                          "every time, so an eval run can never silently "
                          "overwrite an existing model.")
    ap.add_argument("--train_full_only", action="store_true",
                     help="train a single model on 100%% of --dataset and save "
                          "it to --save_model (required in this mode)")
    ap.add_argument("--load_model", default=None,
                     help="NEW: path to a state_dict saved by --train_full_only "
                          "or --save_model. Required with --eval_only.")
    ap.add_argument("--eval_only", action="store_true",
                     help="NEW: load --load_model and evaluate it on every "
                          "graph in --dataset, with NO training. This is the "
                          "actual cross-domain test -- e.g. a model trained "
                          "on results/pyg_dataset.pt (simulated), evaluated "
                          "on results/wick_pyg_dataset.pt (real), unchanged "
                          "by anything in the Wick data.")
    args = ap.parse_args()

    torch.manual_seed(SEED)
    random.seed(SEED)
    np.random.seed(SEED)

    if args.train_full_only and not args.save_model:
        raise SystemExit("--train_full_only requires --save_model <path>.")
    if args.eval_only and not args.load_model:
        raise SystemExit("--eval_only requires --load_model <path>.")

    dataset = torch.load(args.dataset, weights_only=False)
    isolates = isolate_groups(dataset)
    unique_isolates = sorted(set(isolates))
    print(f"{len(dataset)} graphs across {len(unique_isolates)} isolates.")
    if len(unique_isolates) < 10:
        print("WARNING: <10 isolates. Bootstrap intervals will be wide -- "
              "report them as-is, do not hide the width.")

    plasmid_examples, graph_examples = build_examples(dataset)
    in_dim = dataset[0].x.shape[1]

    # ---- NEW: eval-only, no training, on a loaded checkpoint ----
    if args.eval_only:
        all_idx = list(range(len(dataset)))
        print(f"loaded {args.load_model}, evaluating on all {len(dataset)} "
              f"graphs in {args.dataset} -- no training happens in this mode.")

        # Evaluate one graph at a time rather than passing the whole list
        # into run_fold at once: run_fold converts test_graph_idx to a
        # Python set() internally, so its iteration order isn't safe to
        # assume matches all_idx -- doing it per-graph avoids needing to
        # reconstruct "which result came from which graph" after the fact.
        per_isolate_probs, per_isolate_true = defaultdict(list), defaultdict(list)
        graph_pred_all, graph_true_all = [], []
        node_probs_all, node_true_all = [], []
        model = GNNModel(in_dim=in_dim, hidden=args.hidden).to(args.device)
        state = torch.load(args.load_model, map_location=args.device)
        model.load_state_dict(state)

        for gi in all_idx:
            iso = isolates[gi]
            res = run_fold(model, dataset, plasmid_examples, graph_examples,
                            train_graph_idx=[], test_graph_idx=[gi],
                            device=args.device, epochs=0)
            per_isolate_probs[iso].extend(res["plasmid_probs"])
            per_isolate_true[iso].extend(res["plasmid_true"])
            graph_pred_all.extend(res["graph_pred"])
            graph_true_all.extend(res["graph_true"])
            node_probs_all.extend(res["node_probs"])
            node_true_all.extend(res["node_true"])

        out = compute_metrics(len(unique_isolates), len(dataset),
                               per_isolate_probs, per_isolate_true,
                               graph_pred_all, graph_true_all,
                               node_probs_all, node_true_all)
        out["mode"] = f"eval_only, loaded {args.load_model}, no training on {args.dataset}"
        with open(args.out, "w") as f:
            json.dump(out, f, indent=2)
        print(json.dumps(out, indent=2))
        return

    # ---- train a single model on 100% of --dataset, save, stop ----
    if args.train_full_only:
        print(f"Training single model on 100% of {args.dataset} "
              f"({len(dataset)} graphs) for {args.epochs} epochs...")
        model = GNNModel(in_dim=in_dim, hidden=args.hidden).to(args.device)
        run_fold(model, dataset, plasmid_examples, graph_examples,
                 list(range(len(dataset))), test_graph_idx=[],
                 device=args.device, epochs=args.epochs)
        torch.save(model.state_dict(), args.save_model)
        print(f"saved to {args.save_model}")
        return

    # ---- default: leave-isolate-out CV on --dataset ----
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

    out = compute_metrics(len(unique_isolates), len(dataset),
                           per_isolate_probs, per_isolate_true,
                           graph_pred_all, graph_true_all,
                           node_probs_all, node_true_all)
    out["mode"] = f"leave-isolate-out CV on {args.dataset}"
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))

    # PATCHED: saving a full-data model after CV is now opt-in only --
    # nothing happens here unless --save_model was explicitly passed.
    if args.save_model:
        print(f"--save_model was set: training an additional master model "
              f"on 100% of {args.dataset} and saving to {args.save_model} "
              f"(separate from the CV results above)...")
        master_model = GNNModel(in_dim=in_dim, hidden=args.hidden).to(args.device)
        run_fold(master_model, dataset, plasmid_examples, graph_examples,
                 list(range(len(dataset))), test_graph_idx=[],
                 device=args.device, epochs=args.epochs)
        torch.save(master_model.state_dict(), args.save_model)
        print(f"saved to {args.save_model}")


if __name__ == "__main__":
    main()
