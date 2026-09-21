"""
14_visualize_results.py -- turn results/*.json into one readable report:
does the pipeline work, and is any of the headline number trustworthy?

Reads (no torch needed -- this only needs the JSON/CSV side-products):
  results/baselines.json          B0 / B0+ / B1 / B2 (06_baselines.py)
  results/graph_diagnostics.json  per-graph topology (segments/edges/etc.)
  results/gnn_results.json        B3/B4 GraphSAGE (13_train_gnn.py)

Produces results/report.png (4 panels) and prints a plain-text verdict.
Safe to re-run any time results/gnn_results.json is refreshed -- panels
about the confirmatory node/graph heads render as "not computed yet" until
that file has pr_auc_per_class / mae_baseline_predict_mean in it (i.e. until
the patched 13_train_gnn.py has been re-run).
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")


def load(name):
    p = os.path.join(RES, name)
    if not os.path.exists(p):
        return None
    return json.load(open(p))


def panel_baselines(ax, baselines, gnn):
    """Panel A: PR-AUC across every method tried so far, one axis, same
    metric, same slice (with_is_empty = the full, real-world task)."""
    if baselines is None:
        ax.text(0.5, 0.5, "results/baselines.json not found", ha="center")
        return
    slice_ = baselines["with_is_empty"]
    methods = ["B0", "B0_pooled", "B0+", "B1", "B2"]
    vals = [slice_["pr_auc"][m] for m in methods]
    labels = list(methods)
    colors = ["#b0b0b0"] * len(methods)

    if gnn and "b4_plasmid_pooled_head" in gnn:
        methods.append("B4-pooled\n(GraphSAGE)")
        vals.append(gnn["b4_plasmid_pooled_head"]["pr_auc"])
        labels.append("B4-pooled\n(GraphSAGE)")
        colors.append("#d62728")  # flagged red -- see caption

    prevalence = slice_["prevalence"]
    x = np.arange(len(vals))
    ax.bar(x, vals, color=colors)
    ax.axhline(prevalence, ls="--", c="k", lw=1)
    ax.text(len(vals) - 0.5, prevalence + 0.02,
            f"prevalence floor ({prevalence:.2f})", fontsize=8, ha="right")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("PR-AUC")
    ax.set_title("A. Plasmid-level PR-AUC, all methods so far\n"
                  "(red = flagged leaky by the project's own docs, not confirmatory)")


def panel_topology(ax, diag):
    """Panel B: the actual reason B3/B4's node/graph heads are expected to
    struggle -- most nodes have no neighbours to pass messages with."""
    if diag is None:
        ax.text(0.5, 0.5, "results/graph_diagnostics.json not found", ha="center")
        return
    n_graphs = len(diag)
    total_seg = sum(g["n_segments"] for g in diag)
    total_isolated = sum(g["isolated_nodes"] for g in diag)
    total_edges = sum(g["n_simple_edges"] for g in diag)
    edgeless = sum(1 for g in diag if g["n_simple_edges"] == 0)

    labels = ["Isolated\n(no neighbours)", "Connected"]
    vals = [total_isolated, total_seg - total_isolated]
    ax.pie(vals, labels=labels, autopct="%1.0f%%", colors=["#d62728", "#2ca02c"])
    ax.set_title(
        f"B. Node connectivity across all {n_graphs} simulated graphs\n"
        f"{total_edges} real edges total | {edgeless}/{n_graphs} graphs "
        f"({edgeless/n_graphs:.0%}) have zero edges"
    )


def panel_node_head(ax, gnn):
    """Panel C: the actual confirmatory test -- does the node head beat
    just guessing class frequency?"""
    if gnn is None or "confirmatory_node_head" not in gnn:
        ax.text(0.5, 0.5, "no gnn_results.json yet", ha="center")
        return
    nh = gnn["confirmatory_node_head"]
    per_class = nh.get("pr_auc_per_class")
    baseline = nh.get("prevalence_baseline_per_class")
    if not per_class:
        ax.text(0.5, 0.5,
                "NOT YET COMPUTED\n\nresults/gnn_results.json was produced "
                "by the pre-patch script -- it only stored n_examples "
                f"({nh.get('n_examples','?')}), never scored them.\n\n"
                "Re-run the patched 13_train_gnn.py to fill this in.",
                ha="center", va="center", fontsize=9, wrap=True)
        ax.set_title("C. Node head: PR-AUC vs. \"just guess the frequency\"")
        ax.axis("off")
        return
    classes = list(per_class.keys())
    model_v = [per_class[c] if per_class[c] is not None else 0 for c in classes]
    base_v = [baseline[c] if baseline and baseline.get(c) is not None else 0
              for c in classes]
    x = np.arange(len(classes))
    w = 0.35
    ax.bar(x - w / 2, model_v, w, label="GraphSAGE node head", color="#1f77b4")
    ax.bar(x + w / 2, base_v, w, label="prevalence baseline", color="#b0b0b0")
    ax.set_xticks(x)
    ax.set_xticklabels(classes, fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7)
    ax.set_title("C. Node head: PR-AUC vs. \"just guess the frequency\"\n"
                 "(bars this close together = not learning past class balance)")


def panel_graph_head(ax, gnn):
    """Panel D: does the graph head beat 'always predict the mean'?"""
    if gnn is None or "confirmatory_graph_head" not in gnn:
        ax.text(0.5, 0.5, "no gnn_results.json yet", ha="center")
        return
    gh = gnn["confirmatory_graph_head"]
    mae = gh.get("mae_missing_count")
    base = gh.get("mae_baseline_predict_mean")
    if base is None:
        ax.text(0.5, 0.5,
                f"NOT YET COMPUTED\n\nmae_missing_count = {mae:.3f} exists, "
                "but with no baseline to compare it to, that number alone "
                "means nothing.\n\nRe-run the patched 13_train_gnn.py to "
                "get mae_baseline_predict_mean.",
                ha="center", va="center", fontsize=9, wrap=True)
        ax.set_title("D. Graph head: MAE vs. \"always predict the mean\"")
        ax.axis("off")
        return
    x = np.arange(2)
    ax.bar(x, [mae, base], color=["#1f77b4", "#b0b0b0"])
    ax.set_xticks(x)
    ax.set_xticklabels(["GraphSAGE\ngraph head", "always predict\ncohort mean"])
    ax.set_ylabel("MAE (missing-plasmid count)")
    ax.set_title("D. Graph head: MAE vs. \"always predict the mean\"\n"
                 "(lower is better; bars this close = not adding signal)")


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--gnn-results", default="gnn_results.json",
                     help="PATCHED: filename under results/ to read B3/B4 "
                          "numbers from -- pass gnn_results_wick.json once "
                          "that exists, so this can render a report for the "
                          "real-data run without overwriting the simulated one.")
    ap.add_argument("--out", default="report.png",
                     help="Output filename under results/.")
    ap.add_argument("--title-suffix", default="",
                     help="e.g. '(simulated)' or '(Wick, real data)'.")
    args = ap.parse_args()

    baselines = load("baselines.json")
    diag = load("graph_diagnostics.json")
    gnn = load(args.gnn_results)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    panel_baselines(axes[0, 0], baselines, gnn)
    panel_topology(axes[0, 1], diag)
    panel_node_head(axes[1, 0], gnn)
    panel_graph_head(axes[1, 1], gnn)
    suffix = f" {args.title_suffix}" if args.title_suffix else ""
    fig.suptitle(f"Plasmid recovery pipeline -- where the evidence actually stands{suffix}",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(RES, args.out)
    fig.savefig(out, dpi=150)
    print("wrote", out)

    # plain-text verdict, so this doesn't require looking at the image
    print()
    print("=" * 70)
    print("VERDICT")
    print("=" * 70)
    if baselines:
        b2 = baselines["with_is_empty"]["pr_auc"]["B2"]
        print(f"Best feature-based baseline (B2, hand-crafted features, no "
              f"graph): PR-AUC {b2:.3f}")
    if gnn and "b4_plasmid_pooled_head" in gnn:
        print(f"GraphSAGE pooled head (B4, descriptive, leak-flagged): "
              f"PR-AUC {gnn['b4_plasmid_pooled_head']['pr_auc']:.3f} "
              f"-- NOT directly comparable to B2 above; don't headline this.")
    if diag:
        n = len(diag)
        iso_frac = sum(g["isolated_nodes"] for g in diag) / sum(g["n_segments"] for g in diag)
        print(f"Topology: {iso_frac:.0%} of nodes across {n} simulated graphs "
              f"have no edges -- GraphSAGE's message passing has nothing to "
              f"pass for most of the dataset.")
    if gnn and gnn.get("confirmatory_node_head", {}).get("pr_auc_per_class"):
        print("Confirmatory node/graph heads: COMPUTED -- see panels C/D.")
    else:
        print("Confirmatory node/graph heads: NOT YET COMPUTED. The only "
              "number you have (0.97) is the one the pipeline's own code "
              "explicitly flags as not confirmatory. Re-run the patched "
              "13_train_gnn.py before claiming the model works.")


if __name__ == "__main__":
    main()
