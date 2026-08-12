"""
Shared label-leak guard.

Two columns must never reach a feature matrix:

  circular      Sec. IV already excludes Flye's circularity flag, because the
                label definition itself invokes circularity -- a plasmid is
                'recovered' when it comes back as its own (circular) contig, so
                the flag is a restatement of the answer.
  completeness  an assembler's own verdict on whether the replicon came back
                whole. That IS the label, under a different name.

`contigs_summary_sup_cleaned.csv` (NEKSUS, 3,917 rows x 6 assemblers x 92
samples) carries both. It is the one deposited artefact that invites direct
feature extraction -- it is a tidy per-contig table with length, depth, GC-free
metadata already joined -- so the temptation is structural, not hypothetical.
The guard therefore lives here and is called BEFORE the file is ever parsed,
rather than being left as a comment in whichever script gets there first.

Both columns are retained as DIAGNOSTICS. They are legitimate for describing the
corpus and for sanity-checking labels; they are illegitimate as model input.

Usage:
    from leak_guard import assert_no_leak_columns, load_neksus_contigs_summary

    assert_no_leak_columns(FEATS_ALL, "plasmid-head feature matrix")
    feats, diags = load_neksus_contigs_summary(path)   # feats can never carry them
"""
import csv

# Column names banned from any feature matrix, in any casing or separator style.
BANNED = ("circular", "completeness")

# Names that merely CONTAIN a banned token are also refused: 'is_circular',
# 'chromosome_circularity', 'circular_frac', 'completeness_pct' are the same
# leak wearing a different label.
def _normalise(name):
    return str(name).strip().lower().replace("-", "_").replace(" ", "_")


def is_leak_column(name):
    n = _normalise(name)
    return any(b in n for b in BANNED)


def assert_no_leak_columns(columns, where="feature matrix"):
    """Hard-fail if any banned column reached `columns`. Never warns -- a leak
    that only warns gets shipped."""
    bad = [c for c in columns if is_leak_column(c)]
    if bad:
        raise AssertionError(
            "LABEL LEAK BLOCKED in %s: %s. The circularity flag and the "
            "assembler's completeness verdict restate the label; Sec. IV "
            "excludes them. Keep them as diagnostics only." % (where, bad))
    return list(columns)


def load_neksus_contigs_summary(path):
    """Parse the NEKSUS contig summary with the leak columns split off.

    Returns (features, diagnostics):
      features    list of dicts with the banned columns REMOVED. Safe to join
                  into a feature matrix.
      diagnostics list of dicts with only the key columns plus the banned ones.
                  Safe to report; never to train on.

    The split happens at the parser, so no caller can accidentally carry the
    columns forward -- there is no code path that returns them inside `features`.
    """
    features, diagnostics = [], []
    with open(path, newline="", encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        cols = rd.fieldnames or []
        leaks = [c for c in cols if is_leak_column(c)]
        keys = [c for c in cols if c in ("sample", "assembler", "contig_name")]
        for r in rd:
            features.append({k: v for k, v in r.items() if not is_leak_column(k)})
            d = {k: r.get(k) for k in keys}
            d.update({k: r.get(k) for k in leaks})
            diagnostics.append(d)
    assert_no_leak_columns(features[0].keys() if features else [],
                           "NEKSUS contigs_summary features")
    return features, diagnostics


if __name__ == "__main__":
    import sys
    p = sys.argv[1] if len(sys.argv) > 1 else \
        r"D:\DNA-Sequencing\Dataset\29584931\contigs_summary_sup_cleaned.csv"
    feats, diags = load_neksus_contigs_summary(p)
    print("rows                 : %d" % len(feats))
    print("feature columns      : %s" % sorted(feats[0].keys()))
    print("diagnostic-only cols : %s" % sorted(diags[0].keys()))
    for bad in ("circular", "completeness", "no_circular_contigs",
                "chromosome_circularity", "is_circular", "completeness_pct"):
        assert is_leak_column(bad), bad
    try:
        assert_no_leak_columns(["length", "depth", "circular"], "self-test")
        print("SELF-TEST FAILED: guard did not fire")
    except AssertionError:
        print("guard fires on a banned column: OK")
