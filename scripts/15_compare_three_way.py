"""
15_compare_three_way.py -- one chart: simulated (in-domain) vs. B (simulated
model, evaluated on Wick, no retraining) vs. C (trained from scratch on the
7 Wick isolates, leave-isolate-out CV). Same honest-baseline style as
14_visualize_results.py's report.png.

Usage (from repo root):
    python scripts/15_compare_three_way.py \
        --simulated results/gnn_results_simulated.json \
        --wick-transfer results/gnn_results_wick_from_sim_model.json \
        --wick-cv results/gnn_results_wick_cv.json \
        --out results/report_three_way.png
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COND_LABELS = ["Simulated\n(in-domain)", "Sim-model -> Wick\n(no retrain, n=7)", "Wick-only\n(from scratch, n=7)"]
COND_COLORS = ["#888888", "#1f77b4", "#d62728"]


def load(p):
    return json.load(open(p))


def grouped_bars(ax, model_vals, base_vals, title):
    x = np.arange(3)
    w = 0.35
    ax.bar(x - w / 2, model_vals, w, label="model", color=COND_COLORS)
    ax.bar(x + w / 2, base_vals, w, label="baseline", color="#cccccc")
    ax.set_xticks(x)
    ax.set_xticklabels(COND_LABELS, fontsize=8)
    ax.legend(fontsize=7)
    ax.set_title(title, fontsize=10)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--simulated", required=True)
    ap.add_argument("--wick-transfer", required=True)
    ap.add_argument("--wick-cv", required=True)
    ap.add_argument("--out", default="results/report_three_way.png")
    args = ap.parse_args()

    sim = load(args.simulated)
    b = load(args.wick_transfer)
    c = load(args.wick_cv)
    runs = [sim, b, c]

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # A. pooled head PR-AUC with CI whiskers
    ax = axes[0, 0]
    vals = [r["b4_plasmid_pooled_head"]["pr_auc"] for r in runs]
    ci_lo = [r["b4_plasmid_pooled_head"]["bootstrap_ci95"][0] for r in runs]
    ci_hi = [r["b4_plasmid_pooled_head"]["bootstrap_ci95"][1] for r in runs]
    err = [[v - lo for v, lo in zip(vals, ci_lo)], [hi - v for v, hi in zip(vals, ci_hi)]]
    x = np.arange(3)
    ax.bar(x, vals, color=COND_COLORS, yerr=err, capsize=6)
    ax.set_xticks(x); ax.set_xticklabels(COND_LABELS, fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("PR-AUC")
    ax.set_title("A. Pooled head PR-AUC, with 95% CI\n(descriptive/leaky per open item B -- wide CIs on n=7 are real, not an error)", fontsize=10)

    # B. node head: fragmented
    grouped_bars(axes[0, 1],
                 [r["confirmatory_node_head"]["pr_auc_per_class"]["fragmented"] for r in runs],
                 [r["confirmatory_node_head"]["prevalence_baseline_per_class"]["fragmented"] for r in runs],
                 "B. Node head -- fragmented\n(bars above baseline = real signal)")

    # C. node head: recovered
    grouped_bars(axes[1, 0],
                 [r["confirmatory_node_head"]["pr_auc_per_class"]["recovered"] for r in runs],
                 [r["confirmatory_node_head"]["prevalence_baseline_per_class"]["recovered"] for r in runs],
                 "C. Node head -- recovered\n(majority class -- baseline is already high)")

    # D. graph head MAE (lower is better -- flip framing but keep bar style)
    ax = axes[1, 1]
    model_mae = [r["confirmatory_graph_head"]["mae_missing_count"] for r in runs]
    base_mae = [r["confirmatory_graph_head"]["mae_baseline_predict_mean"] for r in runs]
    x = np.arange(3)
    w = 0.35
    ax.bar(x - w / 2, model_mae, w, label="model", color=COND_COLORS)
    ax.bar(x + w / 2, base_mae, w, label="baseline (predict mean)", color="#cccccc")
    ax.set_xticks(x); ax.set_xticklabels(COND_LABELS, fontsize=8)
    ax.legend(fontsize=7)
    ax.set_ylabel("MAE (lower is better)")
    ax.set_title("D. Graph head -- missing-plasmid count MAE", fontsize=10)

    fig.suptitle("Simulated vs. real Wick data -- does the model actually transfer?",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print("wrote", args.out)

    print("\nVERDICT:")
    for name, r in zip(["simulated", "B (transfer)", "C (from-scratch)"], runs):
        frag = r["confirmatory_node_head"]["pr_auc_per_class"]["fragmented"]
        frag_base = r["confirmatory_node_head"]["prevalence_baseline_per_class"]["fragmented"]
        mae = r["confirmatory_graph_head"]["mae_missing_count"]
        mae_base = r["confirmatory_graph_head"]["mae_baseline_predict_mean"]
        beats = "beats" if frag > frag_base and mae < mae_base else "does NOT beat"
        print(f"  {name:20s}: fragmented PR-AUC {frag:.3f} (base {frag_base:.3f}), "
              f"graph MAE {mae:.3f} (base {mae_base:.3f}) -- {beats} its baseline on both")


if __name__ == "__main__":
    main()
