"""
Phase 0 / step 6 -- the baseline ladder, up to B2.

Per Sec. V:
  B0   logistic, length only                  (one intrinsic feature)
  B0+  logistic, length + rel. depth + dGC    (PRIMARY NULL, no topology)
  B1   logistic, all features                 (+ hand-engineered local topology)
  B2   gradient boosting, no edges            (same features, nonlinear)

B3/B4 need the GNN and are not attempted here.

Evaluation follows Sec. V exactly:
  * leave-isolate-out CV, GROUPED so no depth ladder straddles folds
  * PR-AUC as the reported metric (not accuracy, not ROC-AUC)
  * paired bootstrap OVER ISOLATES for model comparisons; DeLong is not used
    because it has no valid form for PR-AUC
  * the 50-95% coverage band is held out and reported separately

A degeneracy is reported alongside, not buried. Sec. III defines S(p) = empty
=> the plasmid is absent, so the is_empty indicator that the plasmid head
receives nearly determines the absent class. Every model is therefore run twice,
with and without is_empty, and is_empty-alone is scored as its own baseline.

--------------------------------------------------------------------------
COHORT v2 CHANGES -- three things this script now refuses to let past.

1. ADEQUACY GATE (runs before any analysis). The matched-pair comparison is
   this arm's primary contribution: same plasmid, same isolate, IDENTICAL
   LENGTH, one prep drops it and the other does not. Inside a matched pair B0
   is at chance BY CONSTRUCTION, so anything beating chance is provably using
   something beyond length. v1 had 5 such pairs; that test cannot run at n=5.
   If fewer than MIN_DISCORDANT plasmids change label across prep or across
   depth, this script prints why the cohort is inadequate and EXITS without
   analysing it. An underpowered cohort does not get scored.

2. n IS ISOLATES, NEVER LABEL UNITS. Depth and prep cells are nested inside
   isolate and contribute no independent degrees of freedom. v1 reported 168
   units when the resampling unit was 12 isolates -- a 14x overstatement of the
   evidence. Every printed n now leads with isolates.

3. UNRESOLVED verdicts. A comparison whose bootstrap interval boundary changes
   sign between the is_empty-in and is_empty-out arms is not 'clearing' in one
   and 'failing' in the other -- it is unresolved, and is now labelled that way.
   v1's B1-B0+ was [+0.002,+0.136] and [-0.000,+0.142]: the same effect either
   side of the threshold on an arm that should not matter.

4. B0 IS THE PLASMID'S OWN LENGTH. It previously selected the pooled mean and
   max log length of the SUPPORTING CONTIGS -- the plasmid's own length was
   never in the design matrix at all -- and that vector is identically zero for
   every absent unit, so the rung designated as the no-topology floor already
   encoded the class it is meant to be ignorant of (v1: 30 of 30 absent units).
   Measured on v1: as coded 0.867, true plasmid-length-only 0.815. Every rung is
   a difference FROM this floor, so it was not a constant offset -- it made
   B0+ - B0, B1 - B0+ and B2 - B0+ uninterpretable. The old definition is kept
   as a separate rung `B0_pooled` so the correction is auditable rather than a
   silent renumbering, and the column choice is asserted at import.

The 'tiny' band (<4 kb) is excluded from the interventional analysis and scored
separately. Below ~4 kb the plasmid is shorter than nearly every read and Flye
drops it in all four cells regardless of prep or depth: those units are
constants, not observations, and pooling them inflates every model on the
ladder.
"""
import os, csv, json, collections, math, sys
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import average_precision_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from leak_guard import assert_no_leak_columns

ROOT = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\work"
RESULTS = os.path.join(ROOT, "results")
RNG = np.random.default_rng(20260811)

MIN_DISCORDANT = 20        # adequacy gate; see (1) above

FEATS_INTRINSIC = ["log10_length", "relative_depth", "delta_gc"]
FEATS_TOPO = ["degree", "component_size", "n_components"]
FEATS_ALL = FEATS_INTRINSIC + FEATS_TOPO + ["depth_x_length"]

# THE PLASMID'S OWN LENGTH -- not a node feature, not pooled.
#
# Every name in FEATS_ALL is a property of a SUPPORTING CONTIG, pooled over S(p)
# as mean and max. `log10_length` there is the log length of the contigs that
# support the plasmid; the plasmid's own length was never in the design matrix.
# So `cols_for(["log10_length"])` -- the old B0 -- was "mean and max log length
# of the supporting contigs", which is (a) not a length-only baseline and (b)
# identically zero for every absent unit, because an absent plasmid has
# S(p) = empty and therefore an all-zero pooled vector. The rung designated as
# the no-topology floor already encoded the class it is meant to be ignorant of.
#
# This sentinel selects a single intrinsic column appended AFTER the pooled
# block and after is_empty, read from `plasmid_len` in plasmid_labels.csv, which
# 04_label.py takes as len(seq) directly from sim/<isolate>_ref.fasta.
#
# It is deliberately NOT added to FEATS_ALL: putting it there would pool it,
# which is the exact bug being fixed.
PLASMID_LEN = "plasmid_log10_length"

# Sec. IV excludes the circularity flag; the same leak reaches this pipeline
# through the NEKSUS contig summary's `circular` / `completeness` columns.
# Assert at the point the feature list is defined, not after it is used.
assert_no_leak_columns(FEATS_ALL, "plasmid-head feature matrix (06_baselines)")
assert_no_leak_columns([PLASMID_LEN], "intrinsic plasmid column (06_baselines)")

FAIL_LABELS = ("absent", "fragmented", "absorbed")
INTERVENTIONAL = "interventional"


# ------------------------------------------------------------------ loading
def load_truth():
    """(isolate, plasmid_name) -> ground-truth row, for stratum / analysis set."""
    p = os.path.join(ROOT, "sim", "ground_truth.csv")
    t = {}
    for r in csv.DictReader(open(p, encoding="utf-8")):
        t[(r["isolate"], r["plasmid_name"])] = r
    return t


def analysis_set_of(gt_row, plasmid_len):
    """Prefer the column written by 02_build_sim_refs.py; derive it from the
    v2 bounds if an older ground_truth.csv is on disk, so this script never
    silently analyses a tiny unit as interventional."""
    if gt_row and gt_row.get("analysis_set"):
        return gt_row["analysis_set"]
    return "deterministic_tiny" if plasmid_len < 4_000 else INTERVENTIONAL


def load():
    nf = collections.defaultdict(dict)
    with open(os.path.join(RESULTS, "node_features.csv")) as fh:
        for r in csv.DictReader(fh):
            nf[r["tag"]][r["segment"]] = {k: float(r[k]) for k in FEATS_ALL}

    truth = load_truth()
    rows = []
    len_mismatch = []
    gdir = os.path.join(ROOT, "graphs")
    for tag in sorted(os.listdir(gdir)):
        lp = os.path.join(gdir, tag, "plasmid_labels.csv")
        if not os.path.exists(lp):
            continue
        # Only the 2x2 cohort (prep x depth, fragment length fixed at
        # 15000,13000). The sim01_d50 pilot used a different fragment length and
        # would confound the prep axis, so it stays in the corpus as a pilot but
        # out of the ladder.
        if "_lig_d" not in tag and "_rap_d" not in tag:
            continue
        parts = tag.split("_")
        iso = parts[0]                       # sim01_lig_d30 -> sim01
        prep = parts[1]
        depth = tag.split("_d")[-1]
        for r in csv.DictReader(open(lp)):
            sup = [s for s in r["support"].split(";") if s]
            g = truth.get((iso, r["plasmid"]))
            L = int(r["plasmid_len"])
            # B0 is now this number and nothing else, so it has to be the right
            # number. plasmid_labels.csv gets it from the reference FASTA and
            # ground_truth.csv from the builder; if they disagree, the labels on
            # disk were produced against a different cohort than sim/ holds.
            if g and int(g["plasmid_len"]) != L:
                len_mismatch.append((iso, r["plasmid"], L, int(g["plasmid_len"])))
            rows.append({
                "tag": tag, "isolate": iso, "prep": prep, "depth": depth,
                "plasmid": r["plasmid"], "plasmid_len": L, "label": r["label"],
                "support": sup,
                "stratum": (g or {}).get("stratum", ""),
                "analysis_set": analysis_set_of(g, L),
                "cointegrated": (g or {}).get("cointegrated", ""),
            })
    if len_mismatch:
        raise SystemExit(
            "STALE LABELS -- plasmid_len disagrees between graphs/*/plasmid_labels.csv\n"
            "and sim/ground_truth.csv for %d unit(s), e.g. %r.\n"
            "The labels on disk were produced against a different cohort. Re-run\n"
            "04_label.py before scoring; do not score this."
            % (len(len_mismatch), len_mismatch[:3]))
    return nf, rows


# ------------------------------------------------------------ adequacy gate
def discordance(rows):
    """Plasmids whose label differs across prep, or across depth, at otherwise
    identical settings.

    Counted over INTERVENTIONAL units only: a tiny plasmid that drops in all
    four cells is a constant, and counting constants toward an adequacy gate
    would defeat the point of the gate.

    Two numbers are produced and both are printed, because they answer
    different questions:
      class-level  the label itself differs (recovered vs fragmented counts)
      binary       recovered-vs-failed flips -- the matched pair the design is
                   actually built around ("one prep drops it")
    The gate is on the class-level count, which is the literal criterion and the
    more permissive of the two, so the binary number is always reported next to
    it rather than hidden.
    """
    idx = {}
    for r in rows:
        if r["analysis_set"] != INTERVENTIONAL:
            continue
        idx[(r["isolate"], r["plasmid"], r["prep"], r["depth"])] = r

    preps = sorted({k[2] for k in idx})
    depths = sorted({k[3] for k in idx}, key=lambda d: int(d))
    units = sorted({(k[0], k[1]) for k in idx})

    def failed(r):
        return r["label"] in FAIL_LABELS

    prep_cls, prep_bin, depth_cls, depth_bin = set(), set(), set(), set()
    pairs_prep = pairs_depth = 0

    for iso, p in units:
        for d in depths:                              # prep contrast, depth held
            got = [idx.get((iso, p, pr, d)) for pr in preps]
            got = [g for g in got if g]
            if len(got) < 2:
                continue
            pairs_prep += 1
            if len({g["label"] for g in got}) > 1:
                prep_cls.add((iso, p))
            if len({failed(g) for g in got}) > 1:
                prep_bin.add((iso, p))
        for pr in preps:                              # depth contrast, prep held
            got = [idx.get((iso, p, pr, d)) for d in depths]
            got = [g for g in got if g]
            if len(got) < 2:
                continue
            pairs_depth += 1
            if len({g["label"] for g in got}) > 1:
                depth_cls.add((iso, p))
            if len({failed(g) for g in got}) > 1:
                depth_bin.add((iso, p))

    return {
        "n_interventional_units": len(units),
        "n_prep_comparisons": pairs_prep,
        "n_depth_comparisons": pairs_depth,
        "prep_discordant_class": len(prep_cls),
        "prep_discordant_binary": len(prep_bin),
        "depth_discordant_class": len(depth_cls),
        "depth_discordant_binary": len(depth_bin),
        "discordant_class": len(prep_cls | depth_cls),
        "discordant_binary": len(prep_bin | depth_bin),
        "isolates_with_discordance": len({i for i, _ in (prep_cls | depth_cls)}),
    }


def adequacy_gate(rows):
    d = discordance(rows)
    print("=" * 68)
    print("ADEQUACY GATE -- discordant matched pairs")
    print("=" * 68)
    print("  interventional plasmid units (tiny excluded) : %d"
          % d["n_interventional_units"])
    print("  prep comparisons made / discordant           : %d / %d  (binary %d)"
          % (d["n_prep_comparisons"], d["prep_discordant_class"],
             d["prep_discordant_binary"]))
    print("  depth comparisons made / discordant          : %d / %d  (binary %d)"
          % (d["n_depth_comparisons"], d["depth_discordant_class"],
             d["depth_discordant_binary"]))
    print("  DISCORDANT PLASMIDS (prep or depth)          : %d   floor %d"
          % (d["discordant_class"], MIN_DISCORDANT))
    print("  of which recovered<->failed flips            : %d"
          % d["discordant_binary"])
    print("  isolates contributing discordance            : %d"
          % d["isolates_with_discordance"])

    if d["discordant_class"] < MIN_DISCORDANT:
        print()
        print("  !! COHORT INADEQUATE -- NOT ANALYSED.")
        print("     %d discordant plasmids is below the floor of %d."
              % (d["discordant_class"], MIN_DISCORDANT))
        print("     The matched-pair comparison is this arm's primary")
        print("     contribution: within a pair the plasmid, the isolate and the")
        print("     LENGTH are identical, so B0 is at chance by construction and")
        print("     anything above chance is provably using something beyond")
        print("     length. That test cannot run at this count.")
        print()
        print("     Rebuild the cohort before scoring it. The knobs that move")
        print("     this number are in 02_build_sim_refs.py (strata bounds,")
        print("     N_ISOLATES, repeat structure) and 07_batch_reads.sh")
        print("     (FRAGLEN -- prep depletion operates around the fragment")
        print("     length, so the band it acts on moves with it).")
        # Written to its own file. A failing gate must not clobber the previous
        # run's baselines.json -- refusing to analyse is not a result, and
        # destroying the last real one on the way out is a second failure.
        gp = os.path.join(RESULTS, "baselines_gate_failed.json")
        json.dump({"gate": "FAILED", "min_discordant": MIN_DISCORDANT,
                   "discordance": d}, open(gp, "w"), indent=2)
        print()
        print("  -> %s  (baselines.json left untouched)" % gp)
        return d, False

    print("  -> PASS")
    return d, True


# ------------------------------------------------------------------ features
def pool(nf_tag, support):
    """mean and max concatenated -- S(p) is an unordered set."""
    if not support:
        return [0.0] * (2 * len(FEATS_ALL)), 1.0
    M = np.array([[nf_tag[s][f] for f in FEATS_ALL] for s in support
                  if s in nf_tag], dtype=float)
    if M.size == 0:
        return [0.0] * (2 * len(FEATS_ALL)), 1.0
    return list(M.mean(axis=0)) + list(M.max(axis=0)), 0.0


def build_matrix(nf, rows, drop_absent=False, subset=INTERVENTIONAL):
    """drop_absent=True gives the clean confirmatory subset.

    Dropping the is_empty COLUMN does not remove the degeneracy: an absent
    plasmid has S(p) = empty, so its whole pooled feature vector is zeros, and a
    model can read absence off the zero vector just as easily as off the
    indicator. The only way to score the comparison without absence carrying it
    is to evaluate on fragmented/absorbed alone -- which is exactly what Sec. V
    already prescribes for the node head ("PR-AUC over fragmented and absorbed
    classes only, since absence has no node").

    subset selects the analysis set: 'interventional' is the cohort proper,
    'deterministic_tiny' is the <4 kb band reported apart, None takes both.
    """
    X, y, groups, meta, gray, dropped = [], [], [], [], 0, 0
    for r in rows:
        if subset is not None and r["analysis_set"] != subset:
            continue
        if r["label"] == "gray_50_95":
            gray += 1
            continue
        if drop_absent and r["label"] == "absent":
            dropped += 1
            continue
        vec, is_empty = pool(nf[r["tag"]], r["support"])
        # column order: [pooled mean | pooled max | is_empty | plasmid length]
        # The last column is the plasmid's OWN length and is never pooled --
        # see PLASMID_LEN above. It is finite for every unit including absent
        # ones, which is the whole point: it does not encode S(p) = empty.
        X.append(vec + [is_empty, math.log10(max(r["plasmid_len"], 1))])
        y.append(1 if r["label"] in FAIL_LABELS else 0)
        groups.append(r["isolate"])
        meta.append(r)
    X = np.array(X)
    if len(X):
        assert X.shape[1] == N_COLS, \
            "design matrix is %d wide, expected %d" % (X.shape[1], N_COLS)
    return X, np.array(y), np.array(groups), meta, gray, dropped


IDX = {f: (i, i + len(FEATS_ALL)) for i, f in enumerate(FEATS_ALL)}
POOLED_COLS = frozenset(range(2 * len(FEATS_ALL)))   # the whole S(p) block
EMPTY_COL = 2 * len(FEATS_ALL)
PLEN_COL = 2 * len(FEATS_ALL) + 1
N_COLS = PLEN_COL + 1


def cols_for(names, with_empty):
    c = []
    for n in names:
        if n == PLASMID_LEN:          # intrinsic: one column, outside the pool
            c.append(PLEN_COL)
            continue
        a, b = IDX[n]                 # pooled: mean and max
        c += [a, b]
    if with_empty:
        c.append(EMPTY_COL)
    return c


def make_model(kind):
    if kind == "gb":
        return HistGradientBoostingClassifier(
            max_depth=3, max_iter=150, learning_rate=0.1,
            min_samples_leaf=5, random_state=0)
    return make_pipeline(StandardScaler(),
                         LogisticRegression(max_iter=5000, C=1.0,
                                            class_weight="balanced"))


def loio_scores(X, y, groups, cols, kind):
    """Leave-isolate-out CV -> out-of-fold score per sample."""
    oof = np.full(len(y), np.nan)
    for iso in np.unique(groups):
        te = groups == iso
        tr = ~te
        if len(np.unique(y[tr])) < 2:
            continue
        m = make_model(kind)
        m.fit(X[tr][:, cols], y[tr])
        oof[te] = m.predict_proba(X[te][:, cols])[:, 1]
    return oof


def prauc(y, s):
    ok = ~np.isnan(s)
    if ok.sum() == 0 or len(np.unique(y[ok])) < 2:
        return float("nan")
    return average_precision_score(y[ok], s[ok])


def paired_bootstrap(y, groups, sa, sb, n=2000):
    """Resample ISOLATES with replacement; recompute both PR-AUCs each time.

    The resampling unit is the isolate, not the label unit: the four
    prep x depth cells of one plasmid are the same plasmid measured four times.
    """
    isos = np.unique(groups)
    diffs = []
    for _ in range(n):
        pick = RNG.choice(isos, size=len(isos), replace=True)
        idx = np.concatenate([np.where(groups == i)[0] for i in pick])
        yy = y[idx]
        if len(np.unique(yy)) < 2:
            continue
        d = prauc(yy, sa[idx]) - prauc(yy, sb[idx])
        if not math.isnan(d):
            diffs.append(d)
    if not diffs:
        return None
    diffs = np.array(diffs)
    return {"mean_diff": float(diffs.mean()),
            "ci_lo": float(np.percentile(diffs, 2.5)),
            "ci_hi": float(np.percentile(diffs, 97.5)),
            "p_gt_0": float((diffs > 0).mean()), "n_boot": len(diffs)}


LADDER = [
    # B0 is the plasmid's own length and NOTHING else -- one column, outside
    # the pooled block. This is the true intrinsic floor every rung is measured
    # from.
    ("B0",        [PLASMID_LEN],    "lr"),
    # The former B0, retained as its own rung rather than deleted: B0_pooled - B0
    # is a direct measurement of how much the contaminated floor was worth, and
    # keeping it is what makes the correction auditable instead of a silent
    # renumbering of every previously reported difference.
    ("B0_pooled", ["log10_length"], "lr"),
    ("B0+",       FEATS_INTRINSIC,  "lr"),
    ("B1",        FEATS_ALL,        "lr"),
    ("B2",        FEATS_ALL,        "gb"),
]
COMPARISONS = [("B1", "B0+"), ("B2", "B0+"), ("B0+", "B0"), ("B0_pooled", "B0")]

# ---------------------------------------------------------------- assertions
# A silent regression here is invisible in the output -- the ladder still prints
# five numbers and they still look plausible. That is how the original bug
# survived an entire cohort. Fail at import, before anything is fitted.
assert PLASMID_LEN not in FEATS_ALL, \
    ("%s must not be in FEATS_ALL -- everything in that list is pooled over "
     "S(p) as mean and max, which is the exact bug this column fixes" % PLASMID_LEN)
_B0_FEATS = next(f for nm, f, _ in LADDER if nm == "B0")
_B0_COLS = cols_for(_B0_FEATS, with_empty=False)
assert len(_B0_COLS) == 1, \
    "B0 must select exactly ONE column (the plasmid's own length); got %r" % (_B0_COLS,)
assert _B0_COLS[0] == PLEN_COL, \
    "B0 must select the intrinsic length column %d; got %d" % (PLEN_COL, _B0_COLS[0])
assert not (set(_B0_COLS) & POOLED_COLS), \
    "B0 column %r falls inside the pooled S(p) block %r" % (_B0_COLS, sorted(POOLED_COLS))
assert PLEN_COL != EMPTY_COL, "intrinsic length column collides with is_empty"
assert set(cols_for(["log10_length"], False)) <= POOLED_COLS, \
    "B0_pooled must be the OLD, pooled definition -- it is the control here"


def run_arm(name, X, y, groups, with_empty, note=""):
    print()
    print("=" * 68)
    print("%s -- leave-isolate-out, PR-AUC" % name)
    if note:
        print(note)
    print("=" * 68)
    n_iso = len(np.unique(groups)) if len(groups) else 0
    base = y.mean() if len(y) else 0.0
    # n IS ISOLATES. The label units below are nested inside these isolates
    # (4 prep x depth cells per plasmid) and carry no independent d.f.
    print("  n = %d ISOLATES  (the resampling and CV unit)" % n_iso)
    print("  nested label units=%d  positives=%d  prevalence=%.3f "
          "(random-classifier PR-AUC)" % (len(y), int(y.sum()), base))
    if len(y) == 0 or y.sum() < 3 or (len(y) - y.sum()) < 3:
        print("  !! too few of one class to fit -- skipped")
        return None

    scores, aucs = {}, {}
    for nm, feats, kind in LADDER:
        s = loio_scores(X, y, groups, cols_for(feats, with_empty), kind)
        scores[nm] = s
        aucs[nm] = prauc(y, s)
        print("  %-9s %-46s %.3f" % (nm, ",".join(feats) +
                                     (" +is_empty" if with_empty else ""), aucs[nm]))
    if with_empty:
        se = loio_scores(X, y, groups, [EMPTY_COL], "lr")
        aucs["is_empty_alone"] = prauc(y, se)
        print("  %-9s %-46s %.3f" % ("DIAG", "is_empty ALONE (degeneracy check)",
                                     aucs["is_empty_alone"]))

    print()
    print("  paired bootstrap over %d isolates (95%% CI on PR-AUC difference):"
          % n_iso)
    boots = {}
    for a, b in COMPARISONS:
        r = paired_bootstrap(y, groups, scores[a], scores[b])
        boots["%s_vs_%s" % (a, b)] = r
        if r:
            verdict = ("EXCEEDS interval" if r["ci_lo"] > 0
                       else "does NOT exceed interval")
            print("    %-13s %+.3f  [%+.3f, %+.3f]  %s"
                  % ("%s-%s" % (a, b), r["mean_diff"], r["ci_lo"], r["ci_hi"], verdict))
    print("    B0_pooled-B0 is the contamination itself, not a model comparison:")
    print("    it is what the old floor bought from pooling S(p) alone. Any value")
    print("    above 0 means every previously reported rung was measured from a")
    print("    floor that already knew the absent class.")
    return {"pr_auc": aucs, "bootstrap": boots, "prevalence": float(base),
            "n_isolates": int(n_iso), "n_units_nested": int(len(y)),
            "n_positive": int(y.sum())}


def cross_arm_stability(arm_a, arm_b):
    """A comparison is UNRESOLVED when a bootstrap interval boundary changes
    sign between the is_empty-in and is_empty-out arms.

    is_empty should not decide the answer: dropping the column leaves an absent
    plasmid with an all-zero pooled vector anyway. So if the verdict moves when
    the column moves, the verdict is a property of the column, not of the
    models -- report it as unresolved rather than as clearing or failing.
    """
    print()
    print("=" * 68)
    print("CROSS-ARM STABILITY (is_empty in vs out)")
    print("=" * 68)
    if not (arm_a and arm_b):
        print("  one arm did not run -- nothing to compare")
        return {}

    def sgn(x):
        return 1 if x > 0 else (-1 if x < 0 else 0)

    out = {}
    for a, b in COMPARISONS:
        k = "%s_vs_%s" % (a, b)
        ra, rb = arm_a["bootstrap"].get(k), arm_b["bootstrap"].get(k)
        if not (ra and rb):
            continue
        moved = [nm for nm in ("ci_lo", "ci_hi") if sgn(ra[nm]) != sgn(rb[nm])]
        clears_a, clears_b = ra["ci_lo"] > 0, rb["ci_lo"] > 0
        if moved or clears_a != clears_b:
            verdict = "UNRESOLVED"
        elif clears_a:
            verdict = "clears in both arms"
        else:
            verdict = "fails in both arms"
        out[k] = {"verdict": verdict, "boundaries_that_changed_sign": moved,
                  "arm_a_ci": [ra["ci_lo"], ra["ci_hi"]],
                  "arm_b_ci": [rb["ci_lo"], rb["ci_hi"]]}
        print("  %-13s A[%+.3f,%+.3f]  B[%+.3f,%+.3f]   %s%s"
              % ("%s-%s" % (a, b), ra["ci_lo"], ra["ci_hi"], rb["ci_lo"],
                 rb["ci_hi"], verdict,
                 ("  (sign change on %s)" % ",".join(moved)) if moved else ""))
    if any(v["verdict"] == "UNRESOLVED" for v in out.values()):
        print()
        print("  UNRESOLVED means exactly that: not 'clears', not 'fails'.")
        print("  Do not report an unresolved comparison as evidence either way.")
    return out


def main():
    nf, rows = load()

    print("=" * 68)
    print("PLASMID-LEVEL DATASET")
    print("=" * 68)
    isolates = sorted({r["isolate"] for r in rows})
    by_set = collections.Counter(r["analysis_set"] for r in rows)
    print("n = %d ISOLATES  (resampling unit; everything below is nested in it)"
          % len(isolates))
    print("nested label units :", len(rows))
    print("analysis sets      :", dict(by_set))
    print("realised classes   :", dict(collections.Counter(r["label"] for r in rows)))
    print("by stratum         :",
          dict(collections.Counter(r["stratum"] for r in rows if r["stratum"])))

    disc, ok = adequacy_gate(rows)
    if not ok:
        return

    X, y, groups, meta, gray, _ = build_matrix(nf, rows)
    print()
    print("interventional subset: %d isolates, %d nested units, %d positives "
          "(%.1f%%)" % (len(np.unique(groups)), len(y), y.sum(),
                        100 * y.mean() if len(y) else 0))
    print("held-out 50-95%% band : %d (reported separately, excluded from fits)"
          % gray)
    if y.sum() < 3 or (len(y) - y.sum()) < 3:
        print("\n!! too few of one class to fit -- stopping")
        return

    out = {"gate": "PASSED", "discordance": disc}
    out["with_is_empty"] = run_arm(
        "ARM A: all mechanisms, is_empty INCLUDED", X, y, groups, True)
    out["without_is_empty"] = run_arm(
        "ARM B: all mechanisms, is_empty EXCLUDED", X, y, groups, False,
        note="  (note: an absent plasmid still has an all-zero pooled vector,\n"
             "   so dropping the column does not remove the degeneracy)")

    # ARM C -- the clean confirmatory subset: absence removed entirely
    Xc, yc, gc, _, grayc, dropped = build_matrix(nf, rows, drop_absent=True)
    out["fragmented_absorbed_only"] = run_arm(
        "ARM C: fragmented/absorbed ONLY (absent dropped)", Xc, yc, gc, False,
        note="  The only arm where absence cannot carry the comparison.\n"
             "  %d absent units removed." % dropped)

    # ARM T -- the <4 kb band, reported apart and never pooled with the above
    Xt, yt, gt, _, grayt, _ = build_matrix(nf, rows, subset="deterministic_tiny")
    out["deterministic_tiny"] = run_arm(
        "ARM T: <4 kb band (REPORTED APART, not an interventional result)",
        Xt, yt, gt, False,
        note="  Below ~4 kb the plasmid is shorter than nearly every read, so\n"
             "  Flye drops it in all four cells regardless of prep or depth.\n"
             "  These units are constants; they are scored for completeness and\n"
             "  are excluded from every claim about prep, depth or topology.")

    out["cross_arm_stability"] = cross_arm_stability(
        out["with_is_empty"], out["without_is_empty"])

    out["dataset"] = {
        "n_isolates": len(isolates),
        "n_units_nested_interventional": int(len(y)),
        "n_positive_interventional": int(y.sum()),
        "gray_band_held_out": gray,
        "analysis_sets": dict(by_set),
        "realised_classes": dict(collections.Counter(r["label"] for r in rows)),
        "realised_classes_interventional": dict(collections.Counter(
            r["label"] for r in rows if r["analysis_set"] == INTERVENTIONAL)),
    }
    json.dump(out, open(os.path.join(RESULTS, "baselines.json"), "w"), indent=2)
    print("\n-> %s" % os.path.join(RESULTS, "baselines.json"))
    print("\nNOTE: n is %d ISOLATES. The %d label units are 4 prep x depth cells"
          % (len(isolates), len(rows)))
    print("per plasmid nested inside those isolates and carry no independent")
    print("degrees of freedom; do not quote them as a sample size.")
    print("\nNOTE: B0 is now the plasmid's own log10 length ALONE -- one column,")
    print("outside the pooled S(p) block. B0+ is still the pooled intrinsic")
    print("triple, so B0+ - B0 spans TWO changes at once: it adds depth and GC,")
    print("and it moves from an intrinsic column to a pooled vector that is zero")
    print("for every absent unit. Read B0_pooled - B0 first to see how much of")
    print("that gap is pooling rather than the added features.")
    print("\nNOTE: B0+ is the primary null. The hypothesis is only testable once")
    print("B4 (the GNN) exists; B1/B2 measure hand-engineered topology and")
    print("nonlinearity, not the learned relational topology the claim is about.")


if __name__ == "__main__":
    main()
