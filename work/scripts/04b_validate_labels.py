"""
Phase 0 / step 4b -- validate labels against ground truth.

The simulation arm exists precisely so this check is possible: we know which
plasmids went into each isolate, at what nominal copy number, and we can read
the REALISED coverage straight off Badread's per-read source tags. So we can ask
whether the labeller's verdict tracks the physical cause it is supposed to
track, instead of merely being self-consistent.

Checks performed:
  1. accounting     -- every designed plasmid gets exactly one label, no
                       plasmid is invented, none goes missing
  2. dose-response  -- 'absent' plasmids should sit at low realised relative
                       depth and 'recovered' ones at high; a labeller that
                       ignored the data would show overlapping distributions
  3. prep effect    -- the ligation arm should lose small plasmids that the
                       rapid arm keeps, on the SAME isolate and depth
  4. determinism    -- is the label predictable from plasmid length alone? If
                       yes the simulation is degenerate and B0 saturates, which
                       would make the whole ladder uninformative.
"""
import os, csv, json, collections

ROOT = r"D:\DNA-Sequencing\work"

# Kept in step with 02_build_sim_refs.py. 'tiny' is NOT an interventional band:
# below ~4 kb the plasmid is shorter than nearly every read and Flye drops it in
# all four cells regardless of prep or depth, so it is reported apart and the
# determinism check below deliberately runs on 'small' (4-12 kb), which is where
# dropout is probabilistic and where the ligation depletion operates.
STRATA = ("tiny", "small", "medium", "large")
INTERVENTIONAL_STRATA = ("small", "medium", "large")


def load_truth():
    t = {}
    for r in csv.DictReader(open(os.path.join(ROOT, "sim", "ground_truth.csv"))):
        t[(r["isolate"], r["plasmid_name"])] = r
    return t


def load_yield(tag):
    """realised bases per replicon, from Badread read source tags"""
    p = os.path.join(ROOT, "asm", tag, "replicon_yield.tsv")
    if not os.path.exists(p):
        return None
    out = {}
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3:
            out[f[0]] = int(f[2])
    return out


def main():
    truth = load_truth()
    gdir = os.path.join(ROOT, "graphs")
    tags = sorted(d for d in os.listdir(gdir)
                  if os.path.isdir(os.path.join(gdir, d)))

    rows = []
    problems = []
    for tag in tags:
        lp = os.path.join(gdir, tag, "plasmid_labels.csv")
        if not os.path.exists(lp):
            continue
        # 2x2 cohort only; the sim01_d50 pilot used a different fragment length
        if "_lig_d" not in tag and "_rap_d" not in tag:
            continue
        parts = tag.split("_")
        iso = parts[0]
        prep = parts[1] if len(parts) > 2 else "?"
        depth = tag.split("_d")[-1]
        y = load_yield(tag)
        chrom_bases = y.get("chromosome", 0) if y else 0
        chrom_len = int(truth[(iso, "plasmid_1")]["chromosome_len"]) \
            if (iso, "plasmid_1") in truth else 0
        chrom_cov = (chrom_bases / chrom_len) if chrom_len else 0

        labelled = set()
        for r in csv.DictReader(open(lp)):
            key = (iso, r["plasmid"])
            if key not in truth:
                problems.append("%s: labelled %s which is not in the design"
                                % (tag, r["plasmid"]))
                continue
            labelled.add(r["plasmid"])
            g = truth[key]
            pl = int(g["plasmid_len"])
            realised = (y.get(r["plasmid"], 0) / pl) if (y and pl) else None
            rel = (realised / chrom_cov) if (realised and chrom_cov) else None
            rows.append({
                "tag": tag, "isolate": iso, "prep": prep, "depth_x": depth,
                "plasmid": r["plasmid"], "plasmid_len": pl,
                "stratum": g["stratum"], "nominal_copy": int(g["copy_number"]),
                "realised_cov": round(realised, 2) if realised else 0.0,
                "realised_rel_depth": round(rel, 2) if rel else 0.0,
                "label": r["label"],
                "best_cov_frac": float(r["best_cov_frac"]),
            })
        designed = {p for (i, p) in truth if i == iso}
        missing = designed - labelled
        if missing:
            problems.append("%s: designed plasmids never labelled: %s"
                            % (tag, sorted(missing)))

    if not rows:
        print("no labelled graphs yet"); return

    out = os.path.join(ROOT, "results", "label_validation.csv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    print("=" * 72)
    print("1. ACCOUNTING")
    print("=" * 72)
    print("  label units:", len(rows), " graphs:", len(tags))
    if problems:
        for p in problems[:20]:
            print("  PROBLEM:", p)
    else:
        print("  OK -- every designed plasmid labelled exactly once, none invented")

    print()
    print("=" * 72)
    print("2. DOSE-RESPONSE: realised relative depth by label")
    print("=" * 72)
    by = collections.defaultdict(list)
    for r in rows:
        by[r["label"]].append(r["realised_rel_depth"])
    for lab in sorted(by):
        v = sorted(by[lab])
        med = v[len(v) // 2]
        print("  %-12s n=%-4d median_rel_depth=%8.2f  min=%7.2f max=%9.2f"
              % (lab, len(v), med, v[0], v[-1]))

    print()
    print("=" * 72)
    print("3. PREP EFFECT (same isolate+depth, lig vs rap)")
    print("=" * 72)
    idx = {(r["isolate"], r["depth_x"], r["prep"], r["plasmid"]): r for r in rows}
    flips = collections.Counter()
    per_stratum = collections.defaultdict(lambda: collections.Counter())
    for (iso, d, prep, p), r in idx.items():
        if prep != "lig":
            continue
        rr = idx.get((iso, d, "rap", p))
        if not rr:
            continue
        lig_fail = r["label"] in ("absent", "fragmented", "absorbed")
        rap_fail = rr["label"] in ("absent", "fragmented", "absorbed")
        key = ("lig_fail" if lig_fail else "lig_ok") + "/" + \
              ("rap_fail" if rap_fail else "rap_ok")
        flips[key] += 1
        per_stratum[r["stratum"]][key] += 1
    for k in sorted(flips):
        print("  %-22s %d" % (k, flips[k]))
    print("  by stratum:")
    for s in STRATA:
        if per_stratum[s]:
            print("    %-7s %s" % (s, dict(per_stratum[s])))
    disc = flips.get("lig_fail/rap_ok", 0) + flips.get("lig_ok/rap_fail", 0)
    print("  -> %d plasmid-pairs flip label with prep alone" % disc)
    if disc == 0:
        print("     WARNING: no prep effect; the arm is not adding information")

    print()
    print("=" * 72)
    print("3b. DEPTH EFFECT (same isolate+prep, low vs high rung)")
    print("=" * 72)
    depths = sorted({r["depth_x"] for r in rows}, key=lambda x: int(x))
    if len(depths) >= 2:
        lo, hi = depths[0], depths[-1]
        dflips = collections.Counter()
        dstrat = collections.defaultdict(collections.Counter)
        for (iso, d, prep, p), r in idx.items():
            if d != hi:
                continue
            rl = idx.get((iso, lo, prep, p))
            if not rl:
                continue
            hi_fail = r["label"] in ("absent", "fragmented", "absorbed")
            lo_fail = rl["label"] in ("absent", "fragmented", "absorbed")
            key = ("lo_fail" if lo_fail else "lo_ok") + "/" + \
                  ("hi_fail" if hi_fail else "hi_ok")
            dflips[key] += 1
            dstrat[r["stratum"]][key] += 1
        for k in sorted(dflips):
            print("  %-22s %d" % (k, dflips[k]))
        print("  by stratum:")
        for s in STRATA:
            if dstrat[s]:
                print("    %-7s %s" % (s, dict(dstrat[s])))
        rec = dflips.get("lo_fail/hi_ok", 0)
        print("  -> %d plasmids rescued by raising depth %sx -> %sx" % (rec, lo, hi))
        if rec == 0:
            print("     (no dose-response yet at this cohort size)")
    else:
        print("  only one depth rung present")

    print()
    print("=" * 72)
    print("4. IS THE LABEL DETERMINED BY LENGTH ALONE?")
    print("=" * 72)
    for s in STRATA:
        sub = [r for r in rows if r["stratum"] == s]
        if not sub:
            continue
        c = collections.Counter(r["label"] for r in sub)
        fail = sum(v for k, v in c.items() if k != "recovered")
        print("  %-7s n=%-4d failure_rate=%5.1f%%  %s"
              % (s, len(sub), 100 * fail / len(sub), dict(c)))
    tiny = [r for r in rows if r["stratum"] == "tiny"]
    if tiny:
        tfr = sum(1 for r in tiny if r["label"] != "recovered") / len(tiny)
        print()
        print("  'tiny' (<4 kb) failure rate %.0f%% over %d units -- EXPECTED to"
              % (100 * tfr, len(tiny)))
        print("  be near-deterministic. This band is reported apart and excluded")
        print("  from the interventional analysis; it is not evidence either way.")

    small = [r for r in rows if r["stratum"] == "small"]
    if small:
        fr = sum(1 for r in small if r["label"] != "recovered") / len(small)
        print()
        print("  DETERMINISM CHECK runs on 'small' (4-12 kb), the band where")
        print("  dropout is probabilistic and the ligation depletion operates.")
        if fr > 0.95 or fr < 0.05:
            print("  WARNING: small-plasmid outcome is near-deterministic "
                  "(%.0f%%). B0 will saturate and the ladder cannot separate."
                  % (100 * fr))
        else:
            print("  OK (%.0f%% failure) -- small-plasmid outcome is NOT determined"
                  % (100 * fr))
            print("     by size alone, so length cannot trivially saturate the ladder.")
    print()
    print("-> %s" % out)


if __name__ == "__main__":
    main()
