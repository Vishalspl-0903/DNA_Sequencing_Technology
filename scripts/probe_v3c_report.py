"""
v3 DESIGN PROBE ROUND 3 -- read out the length x copy-number map and SET the
stratum from it.

The cell-spread gate is the failure this exists to fix. It is computed as the
range of failure rates across the four (prep, depth) cells over interventional
units, and v2 scored 5.5 pp against a 10 pp floor because its 'small' stratum
(4-12 kb x copy 8-30) mixed three regimes: units that fail in both preps, units
that pass in both, and the thin band that actually flips.

A unit only contributes to the spread if it FLIPS on prep. So this script:

  1. prints the label of every plasmid in all five cells;
  2. finds which (length, copy-number) cells flip on prep at BOTH depths --
     flipping at one depth only is not robust enough to build a stratum on;
  3. reports what fraction of a stratum drawn from that region would flip, and
     therefore what cell spread to expect. That number is the prediction the
     v3 cohort has to beat, stated before the cohort is built.

It deliberately does NOT average over the whole probe: an average over a badly
chosen band is exactly the mistake v2 made.
"""
import os, csv, collections

ROOT  = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\work"
PROBE = os.path.join(ROOT, "probe3")
FAIL = ("absent", "fragmented", "absorbed")

CELLS = [("probe03_lig_d30", "lig", 30), ("probe03_rap_d30", "rap", 30),
         ("probe03_lig_d15", "lig", 15), ("probe03_rap_d15", "rap", 15)]
R9 = ("probe03r9_rap_d15", "rap", 15)


def load(tag):
    p = os.path.join(PROBE, "graphs", tag, "plasmid_labels.csv")
    if not os.path.exists(p):
        return None
    return {r["plasmid"]: r["label"] for r in csv.DictReader(open(p))}


def main():
    truth = list(csv.DictReader(open(os.path.join(PROBE, "sim", "probe_truth.csv"))))
    truth.sort(key=lambda r: (int(r["copy_number"]), int(r["plasmid_len"])))
    lab = {t: load(t) for t, _, _ in CELLS}
    r9lab = load(R9[0])
    missing = [t for t in lab if lab[t] is None]
    if missing:
        raise SystemExit("no labels for %r -- run probe_v3c.sh then probe_v3c_label.sh"
                         % missing)

    def failed(tag, name):
        l = lab[tag].get(name)
        return None if l is None else (l in FAIL)

    print("=" * 92)
    print("LENGTH x COPY NUMBER, ALL FOUR GATE CELLS   (F = failed, . = recovered)")
    print("=" * 92)
    print("  %-8s %6s | %5s %5s | %5s %5s | %-14s %s"
          % ("len", "copy", "L30", "R30", "L15", "R15", "flips on prep", "r9 R15"))
    flip_rows = []
    for t in truth:
        n, L, cp = t["plasmid_name"], int(t["plasmid_len"]), int(t["copy_number"])
        f = {tag: failed(tag, n) for tag, _, _ in CELLS}
        m = lambda v: "?" if v is None else ("F" if v else ".")
        f30 = f["probe03_lig_d30"] != f["probe03_rap_d30"]
        f15 = f["probe03_lig_d15"] != f["probe03_rap_d15"]
        note = ("BOTH depths" if (f30 and f15) else
                "30x only" if f30 else "15x only" if f15 else "")
        r9 = r9lab.get(n, "?") if r9lab else "-"
        r9m = "F" if r9 in FAIL else ("." if r9 != "?" else "?")
        print("  %-8d %6d | %5s %5s | %5s %5s | %-14s %s"
              % (L, cp, m(f["probe03_lig_d30"]), m(f["probe03_rap_d30"]),
                 m(f["probe03_lig_d15"]), m(f["probe03_rap_d15"]), note, r9m))
        flip_rows.append({"len": L, "copy": cp, "both": f30 and f15,
                          "any": f30 or f15,
                          "fail_all": all(v for v in f.values() if v is not None),
                          "pass_all": not any(v for v in f.values() if v is not None)})

    both = [r for r in flip_rows if r["both"]]
    anyf = [r for r in flip_rows if r["any"]]
    print()
    print("  flips at BOTH depths : %d / %d" % (len(both), len(flip_rows)))
    print("  flips at either      : %d / %d" % (len(anyf), len(flip_rows)))
    print("  fails in all 4 cells : %d  (constants)" % sum(1 for r in flip_rows if r["fail_all"]))
    print("  passes in all 4 cells: %d  (constants)" % sum(1 for r in flip_rows if r["pass_all"]))

    print()
    print("=" * 92)
    print("WHERE THE FLIPS LIVE")
    print("=" * 92)
    src = both or anyf
    if not src:
        print("  NO plasmid flips on prep anywhere in %d-%d bp."
              % (min(r["len"] for r in flip_rows), max(r["len"] for r in flip_rows)))
        print("  The interventional stratum cannot be set from this probe.")
        return
    tag = "BOTH depths" if both else "either depth (weaker evidence)"
    lo, hi = min(r["len"] for r in src), max(r["len"] for r in src)
    cps = sorted({r["copy"] for r in src})
    print("  using flips at %s" % tag)
    print("  length : %d - %d bp" % (lo, hi))
    print("  copy   : %s" % cps)
    by_cp = collections.Counter(r["copy"] for r in src)
    print("  flips per copy level: %s" % dict(sorted(by_cp.items())))

    # what a stratum drawn from that box would actually deliver
    box = [r for r in flip_rows if lo <= r["len"] <= hi and r["copy"] in cps]
    rate = len(src) / len(box) if box else 0.0
    print()
    print("  a stratum drawn from that box holds %d probe plasmids, of which"
          % len(box))
    print("  %d flip -> %.0f%% flip rate." % (len(src), 100 * rate))
    print()
    print("=" * 92)
    print("PREDICTED CELL SPREAD FOR v3   (state this BEFORE building the cohort)")
    print("=" * 92)
    for share in (0.40, 0.50, 0.60):
        print("  if %.0f%% of interventional units come from that box, and the box"
              % (100 * share))
        print("     flips at %.0f%%, expected spread ~ %.0f pp   (floor 10 pp)"
              % (100 * rate, 100 * share * rate))
    print()
    print("  v2 for comparison: 5.5 pp. Its 'small' band was 4-12 kb x copy 8-30,")
    print("  in which the flip rate measured here is %.0f%% over the whole band."
          % (100 * len(anyf) / len(flip_rows)))

    print()
    print("=" * 92)
    print("SET IN 02_build_sim_refs.py")
    print("=" * 92)
    print('  STRATUM_BOUNDS["small"] = (%d, %d)' % (lo, hi + 1))
    print("  COPY_RANGE[\"small\"]     = (%d, %d)" % (min(cps), max(cps)))
    print("  and raise the 'small' share of COMPOSITIONS to >= 50%")


if __name__ == "__main__":
    main()
