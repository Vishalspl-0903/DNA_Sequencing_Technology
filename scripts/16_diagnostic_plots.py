"""
16_diagnostic_plots.py -- the standard battery of plots used to sanity-check
a sequence/graph model, applied to this project's GraphSAGE checkpoint:

  1. UMAP (falls back to t-SNE if umap-learn isn't installed) of node
     embeddings, colored by each node's plasmid-support label.
  2. Force-directed layout of each assembly graph (networkx spring_layout),
     colored the same way -- lets you SEE which nodes GraphSAGE actually had
     neighbours to pass messages with, not just read the isolated-node %.
  3. Precision-recall curves per class for the node head (the actual curve,
     not just the AP/PR-AUC scalar already in the JSON results).
  4. Feature importance: permutation importance on the 8 raw node features
     (always computed, no extra dependency) -- how much does the node
     head's PR-AUC drop if this feature is shuffled? Also attempts a real
     SHAP summary plot via shap.KernelExplainer if the `shap` package is
     installed; falls back to permutation-importance-only with a clear
     note if not, rather than failing.

This script loads a trained checkpoint and re-runs inference itself (not
just reading the summary JSON), since the plots need raw per-example
predictions and embeddings that 13_train_gnn.py's JSON output doesn't keep.

Usage (from repo root, same env as 13_train_gnn.py -- torch/torch_geometric
required):
    python scripts/16_diagnostic_plots.py \
        --dataset results/wick_pyg_dataset.pt \
        --model results/gnn_model_simulated.pt \
        --out-dir results/diagnostics_wick \
        --title-suffix "(Wick, sim-trained model)"

Run it once per dataset/model combination you want plots for -- e.g. once
for the simulated in-domain run, once for the Wick transfer run.
"""
import argparse
import importlib.util
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn.functional as F
import networkx as nx
from sklearn.metrics import precision_recall_curve, average_precision_score
from sklearn.manifold import TSNE

FEATS = ["log10_length", "relative_depth", "delta_gc", "degree",
         "component_size", "n_components", "depth_x_length"]
NODE_CLASSES = ["fragmented", "absorbed", "recovered"]


def load_train_module():
    """13_train_gnn.py can't be `import`ed by name (starts with a digit) --
    load it by file path instead so this script reuses the exact same
    GNNModel / build_examples / isolate_groups code the training run used,
    rather than a second, possibly-drifted copy of the model definition."""
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "13_train_gnn.py")
    spec = importlib.util.spec_from_file_location("train_gnn", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_inference(mod, dataset, model, device="cpu"):
    """One eval-mode, no-grad pass over every graph, collecting everything
    the plots need: node embeddings, node-head predictions/true labels
    (per plasmid unit with non-empty support), and which graph/isolate each
    node belongs to (for coloring/labeling)."""
    plasmid_examples, _ = mod.build_examples(dataset)
    model.eval()
    all_node_emb, all_node_color, all_node_graph = [], [], []
    node_probs, node_true = [], []

    with torch.no_grad():
        for gi, data in enumerate(dataset):
            d = data.to(device)
            emb = model.encoder(d.x, d.edge_index)
            # default color: "unlabeled" (not in any plasmid's support set)
            node_label = ["unlabeled"] * d.x.shape[0]
            for gi2, unit in plasmid_examples:
                if gi2 != gi:
                    continue
                for idx in unit["support_nodes"]:
                    node_label[idx] = unit["label"]
                if len(unit["support_nodes"]) > 0 and unit["label"] in NODE_CLASSES:
                    nlogits = model.node_head(emb, unit["support_nodes"])
                    nprob = F.softmax(nlogits, dim=1)[0].cpu().numpy()
                    node_probs.append(nprob)
                    node_true.append(NODE_CLASSES.index(unit["label"]))
            all_node_emb.append(emb.cpu().numpy())
            all_node_color.extend(node_label)
            all_node_graph.extend([dataset[gi].tag if hasattr(dataset[gi], "tag") else str(gi)] * d.x.shape[0])

    return {
        "node_emb": np.concatenate(all_node_emb, axis=0),
        "node_color": all_node_color,
        "node_graph": all_node_graph,
        "node_probs": np.array(node_probs) if node_probs else np.zeros((0, 3)),
        "node_true": np.array(node_true),
    }


LABEL_COLORS = {
    "recovered": "#2ca02c", "fragmented": "#d62728",
    "absorbed": "#ff7f0e", "absent": "#7f7f7f", "unlabeled": "#cccccc",
}


def plot_embedding(res, out_dir, suffix):
    """Panel 1: UMAP (or t-SNE fallback) of node embeddings."""
    emb = res["node_emb"]
    if emb.shape[0] < 5:
        print("skipping embedding plot -- too few nodes"); return
    try:
        import umap
        n_neighbors = min(15, max(2, emb.shape[0] // 5))
        # PATCHED: default n_neighbors=15 is ~20% of a 72-node dataset --
        # far too large a fraction for data this small, and can force
        # spurious global structure (e.g. one cluster wrongly pulled far
        # from everything else). Scale down for small inputs, same
        # principle already applied to the t-SNE fallback's perplexity.
        reducer = umap.UMAP(random_state=20260812, n_neighbors=n_neighbors)
        method = f"UMAP (n_neighbors={n_neighbors})"
    except ImportError:
        reducer = TSNE(n_components=2, random_state=20260812,
                        perplexity=min(30, max(2, emb.shape[0] // 3)))
        method = "t-SNE (umap-learn not installed -- pip install umap-learn for the real thing)"
    xy = reducer.fit_transform(emb)

    fig, ax = plt.subplots(figsize=(8, 7))
    colors = res["node_color"]
    for label in set(colors):
        mask = [c == label for c in colors]
        ax.scatter(xy[mask, 0], xy[mask, 1], s=25,
                   c=LABEL_COLORS.get(label, "#000000"), label=label, alpha=0.8)
    ax.legend(fontsize=8)
    ax.set_title(f"Node embeddings, {method} {suffix}")
    fig.tight_layout()
    p = os.path.join(out_dir, "1_embedding_projection.png")
    fig.savefig(p, dpi=150); plt.close(fig)
    print("wrote", p)


def plot_force_directed(dataset, res, out_dir, suffix, max_graphs=12):
    """Panel 2: force-directed layout of each assembly graph, colored by
    node label -- literally shows what GraphSAGE had to work with."""
    n = min(len(dataset), max_graphs)
    cols = min(4, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
    axes = np.array(axes).reshape(-1)

    offset = 0
    for i, data in enumerate(dataset[:n]):
        ax = axes[i]
        nnodes = data.x.shape[0]
        colors = res["node_color"][offset:offset + nnodes]
        offset += nnodes
        G = nx.Graph()
        G.add_nodes_from(range(nnodes))
        ei = data.edge_index.numpy()
        G.add_edges_from(zip(ei[0], ei[1]))
        pos = nx.spring_layout(G, seed=20260812)
        node_colors = [LABEL_COLORS.get(c, "#000000") for c in colors]
        nx.draw(G, pos, ax=ax, node_color=node_colors, node_size=120,
                with_labels=False, edge_color="#aaaaaa")
        tag = getattr(data, "tag", f"graph {i}")
        ax.set_title(tag, fontsize=8)
    for j in range(n, len(axes)):
        axes[j].axis("off")

    fig.suptitle(f"Assembly graphs, force-directed layout {suffix}", fontsize=12)
    fig.tight_layout()
    p = os.path.join(out_dir, "2_force_directed_graphs.png")
    fig.savefig(p, dpi=150); plt.close(fig)
    print("wrote", p)


def plot_pr_curves(res, out_dir, suffix):
    """Panel 3: real precision-recall curves per class, not just the AP
    scalar -- shows where on the curve the model actually sits."""
    if res["node_probs"].shape[0] == 0:
        print("skipping PR curves -- no node-head examples"); return
    fig, ax = plt.subplots(figsize=(7, 6))
    for ci, cname in enumerate(NODE_CLASSES):
        y = (res["node_true"] == ci).astype(int)
        if len(set(y)) < 2:
            continue
        p, r, _ = precision_recall_curve(y, res["node_probs"][:, ci])
        ap = average_precision_score(y, res["node_probs"][:, ci])
        ax.plot(r, p, label=f"{cname} (AP={ap:.3f})")
    ax.set_xlabel("Recall"); ax.set_ylabel("Precision")
    ax.set_ylim(0, 1.05); ax.set_xlim(0, 1.02)
    ax.legend()
    ax.set_title(f"Node head precision-recall curves {suffix}")
    fig.tight_layout()
    p = os.path.join(out_dir, "3_pr_curves.png")
    fig.savefig(p, dpi=150); plt.close(fig)
    print("wrote", p)


def plot_feature_importance(mod, dataset, model, out_dir, suffix, device="cpu"):
    """Panel 4: permutation importance on the raw node features (always
    computed), plus a real SHAP summary plot if `shap` is installed."""
    in_dim = dataset[0].x.shape[1]
    plasmid_examples, _ = mod.build_examples(dataset)
    units = [(gi, u) for gi, u in plasmid_examples
             if len(u["support_nodes"]) > 0 and u["label"] in NODE_CLASSES]
    if not units:
        print("skipping feature importance -- no node-head examples"); return

    def score_with_permutation(feat_idx, rng):
        y_true, y_prob = [], []
        with torch.no_grad():
            for gi, unit in units:
                d = dataset[gi].clone()
                if feat_idx is not None:
                    perm = rng.permutation(d.x.shape[0])
                    d.x[:, feat_idx] = d.x[perm, feat_idx]
                emb = model.encoder(d.x, d.edge_index)
                nlogits = model.node_head(emb, unit["support_nodes"])
                nprob = F.softmax(nlogits, dim=1)[0].numpy()
                y_true.append(NODE_CLASSES.index(unit["label"]))
                y_prob.append(nprob)
        y_true = np.array(y_true); y_prob = np.array(y_prob)
        aps = [average_precision_score((y_true == c).astype(int), y_prob[:, c])
               for c in range(3) if len(set((y_true == c).astype(int))) > 1]
        return float(np.mean(aps)) if aps else float("nan")

    rng = np.random.default_rng(20260812)
    # PATCHED: a single permutation draw on n=26 examples is noisy enough to
    # show a feature as "harmful" (negative drop) purely by chance. Average
    # over several draws so the bars reflect real importance, not one
    # unlucky/lucky shuffle.
    n_repeats = 8
    baseline = score_with_permutation(None, rng)
    drops = {fname: [] for fname in FEATS}
    for _ in range(n_repeats):
        for fi, fname in enumerate(FEATS):
            shuffled = score_with_permutation(fi, rng)
            drops[fname].append(baseline - shuffled)
    drop_mean = {k: float(np.mean(v)) for k, v in drops.items()}
    drop_std = {k: float(np.std(v)) for k, v in drops.items()}

    fig, ax = plt.subplots(figsize=(8, 5))
    names = sorted(drop_mean, key=lambda k: drop_mean[k])
    ax.barh(names, [drop_mean[n] for n in names],
            xerr=[drop_std[n] for n in names], capsize=3, color="#1f77b4")
    ax.axvline(0, c="k", lw=1)
    ax.set_xlabel(f"Mean PR-AUC drop when shuffled ({n_repeats} repeats, error bars = std)")
    ax.set_title(f"Permutation feature importance {suffix}\n"
                 f"(baseline mean PR-AUC = {baseline:.3f}, n={len(units)} examples -- "
                 f"treat bars smaller than their error bar as noise, not signal)")
    fig.tight_layout()
    p = os.path.join(out_dir, "4_feature_importance.png")
    fig.savefig(p, dpi=150); plt.close(fig)
    print("wrote", p)

    try:
        import shap
        n_feat = len(FEATS)
        X = np.stack([dataset[gi].x[unit["support_nodes"][0]][:n_feat].numpy()
                       for gi, unit in units])  # one representative node per unit

        def predict_fn(X_batch):
            """SHAP needs a raw-features -> prediction function. This
            approximates each row as an isolated single node (no edges) --
            i.e. what the node head would predict from that node's own
            features alone, with no graph context. A real per-graph
            explanation would need the full graph each node sits in, which
            KernelExplainer's perturbation approach can't preserve; this is
            a stated simplification, not a bug. Pads a 0 for the model's
            extra "is-global-node" input column, which real segment nodes
            (as opposed to the graph's virtual global node) always have as
            0 anyway."""
            out = []
            with torch.no_grad():
                for row in X_batch:
                    padded = np.concatenate([row, [0.0] * (in_dim - n_feat)])
                    x = torch.tensor(padded, dtype=torch.float32).unsqueeze(0)
                    edge_index = torch.zeros((2, 0), dtype=torch.long)
                    emb = model.encoder(x, edge_index)
                    nlogits = model.node_head(emb, [0])
                    out.append(F.softmax(nlogits, dim=1).numpy()[0])
            return np.array(out)

        n_bg = min(20, len(X))
        explainer = shap.KernelExplainer(predict_fn, X[:n_bg])
        raw = explainer.shap_values(X[:n_bg], nsamples=100)

        # PATCHED AGAIN: different shap versions return multi-class output
        # differently -- older versions: a list of 3 arrays, each
        # (n_samples, n_features). Newer versions: one ndarray, shape
        # (n_samples, n_features, n_classes). Guessing which one silently
        # produced only one (wrong) plot last time. Detect explicitly and
        # normalize to a plain list of clean 2D arrays before doing anything
        # else, so the rest of this code never has to guess again.
        print(f"shap raw output: type={type(raw)}, "
              f"shape={getattr(raw, 'shape', [getattr(a, 'shape', '?') for a in raw] if isinstance(raw, list) else '?')}")
        if isinstance(raw, list):
            per_class = raw                                   # old convention
        elif isinstance(raw, np.ndarray) and raw.ndim == 3:
            per_class = [raw[:, :, c] for c in range(raw.shape[2])]  # new convention
        elif isinstance(raw, np.ndarray) and raw.ndim == 2:
            per_class = [raw]                                 # single-output, unexpected here but handled
        else:
            per_class = None

        if per_class is None:
            print(f"Unrecognized shap output shape -- skipping SHAP plots. "
                  f"Paste the 'shap raw output' line above back to me.")
        else:
            for ci, cname in enumerate(NODE_CLASSES):
                if ci >= len(per_class):
                    continue
                sv = per_class[ci]
                if sv.ndim != 2 or sv.shape[1] != len(FEATS):
                    print(f"skipping {cname}: unexpected shape {sv.shape}")
                    continue
                shap.summary_plot(sv, X[:n_bg], feature_names=FEATS, show=False)
                fig = plt.gcf()
                p = os.path.join(out_dir, f"4b_shap_summary_{cname}.png")
                fig.savefig(p, dpi=150, bbox_inches="tight")
                plt.close(fig)
                print("wrote", p)

            if len(per_class) > 1:
                shap.summary_plot(per_class, X[:n_bg], feature_names=FEATS,
                                   plot_type="bar", class_names=NODE_CLASSES, show=False)
                fig = plt.gcf()
                p = os.path.join(out_dir, "4c_shap_summary_bar_all_classes.png")
                fig.savefig(p, dpi=150, bbox_inches="tight")
                plt.close(fig)
                print("wrote", p)
    except ImportError:
        print("shap not installed (pip install shap) -- skipped SHAP summary "
              "plot, permutation importance above is the substitute.")
    except Exception as e:
        print(f"SHAP summary plot failed ({e}) -- permutation importance "
              f"above is still valid, this is just a bonus panel.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--title-suffix", default="")
    ap.add_argument("--hidden", type=int, default=64)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    mod = load_train_module()
    dataset = torch.load(args.dataset, weights_only=False)
    in_dim = dataset[0].x.shape[1]
    model = mod.GNNModel(in_dim=in_dim, hidden=args.hidden)
    model.load_state_dict(torch.load(args.model, map_location="cpu"))

    res = run_inference(mod, dataset, model)
    suffix = args.title_suffix

    plot_embedding(res, args.out_dir, suffix)
    plot_force_directed(dataset, res, args.out_dir, suffix)
    plot_pr_curves(res, args.out_dir, suffix)
    plot_feature_importance(mod, dataset, model, args.out_dir, suffix)
    print(f"\nAll diagnostic plots written to {args.out_dir}/")


if __name__ == "__main__":
    main()
