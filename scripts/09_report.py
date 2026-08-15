"""
Phase 0 / step 9 -- assemble everything Phase 0 is required to report.

Sec. VII lists what Phase 0 must report before ANY model comparison is
interpreted:
  realised class distribution across the four label classes; node and edge
  counts per graph; number of isolates falling back to the fragmented
  chromosome rule; self-loop counts; simulated-versus-real graph statistic
  comparison; per-mode Flye recovery from the correctness gate.

Everything available from the simulation arm is produced here. The two items
that require real cohorts (the correctness gate, and simulated-vs-real graph
statistics) are marked PENDING with the reason, not quietly omitted.

COHORT v2 ADDITION -- CLASS-COUNT FLOORS AND THE REBUILD GATE PANEL.

v1 realised fragmented 8 and absorbed 1. A mechanism observed once is not a
mechanism that was tested, so any per-class performance claim about it would be
a claim about a single event. Section 0 below prints every realised class
against its floor and states plainly whether the floor is met; classes below
their floor are marked NOT TESTABLE and no per-class performance claim may be
made for them.

Section 0b prints the four gates that decide whether the cohort proceeds to
B3/B4 at all. All four must pass or the cohort is rebuilt.
"""
import os, csv, json, collections, statistics, tarfile

ROOT = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\work"
RES = os.path.join(ROOT, "results")

# Floors on realised per-class counts. Below the floor, the class exists in the
# corpus but cannot carry a performance claim.
#
# 'fragmented' HAD a floor of 20 and no longer has one. It is not a relaxation
# after a near miss (v2 realised 19); the floor was withdrawn because three
# probes established the class is NOT INDUCIBLE in this simulator, so gating on
# it would gate on something no cohort can deliver. Measured, work/probe*:
#
#   cohort v2      IS carrier -> fragmented   0 / 148   (elements 1.0-1.9 kb)
#   probe round 1  repeat 3, 8, 15, 25 kb at 99 % identity, 2 copies -> 0 / 4
#   probe round 2  repeat 2-25 kb at 100 % identity (verbatim-verified),
#                  2-6 copies in plasmid + 2-4 in chromosome, r10_hq @ 30x
#                  AND r9_raw @ 15x                              -> 0 / 12
#
# Nor is it depth: in v2, medium/large units below 1.5x relative depth came out
# 56 recovered vs 6 fragmented, and fragmented sits at a HIGHER median relative
# depth than recovered (3.99 vs 2.87). Badread samples uniformly from a perfect
# circular template and Flye traverses that circle correctly whatever is inside
# it; real fragmentation comes from coverage dropouts, chimeras and strain
# heterogeneity the simulator does not reproduce.
#
# The class is still counted and printed. It is simply not a rebuild trigger,
# and no per-class claim rests on it -- that claim moves to the real cohorts.
CLASS_FLOORS = {"absorbed": 20}
UNGATED_CLASSES = {
    "fragmented": "not inducible in this simulator (3 probes, 0/164); "
                  "reported only, claim moves to the real cohorts",
}
MIN_DISCORDANT = 20
# Cell failure rates must SPREAD, not sit flat. v1 was 26.2/21.4/21.4/23.8 --
# a 4.8-point range across the four cells, i.e. prep and depth were doing
# nothing at the label level while the 17x depletion measured at the DEPTH level
# never reached it. Anything under this range means the interventional axes are
# not moving the outcome.
MIN_CELL_SPREAD_PP = 10.0

INTERVENTIONAL = "interventional"


def section(t):
    print()
    print("=" * 74)
    print(t)
    print("=" * 74)


def load_truth():
    p = os.path.join(ROOT, "sim", "ground_truth.csv")
    if not os.path.exists(p):
        return {}
    return {(r["isolate"], r["plasmid_name"]): r
            for r in csv.DictReader(open(p, encoding="utf-8"))}


def analysis_set_of(g, plasmid_len):
    if g and g.get("analysis_set"):
        return g["analysis_set"]
    return "deterministic_tiny" if plasmid_len < 4_000 else INTERVENTIONAL


def collect_labels(truth):
    """Every labelled plasmid unit in the 2x2 cohort, with its analysis set."""
    rows = []
    gdir = os.path.join(ROOT, "graphs")
    for t in sorted(d for d in os.listdir(gdir)
                    if os.path.isdir(os.path.join(gdir, d))):
        p = os.path.join(gdir, t, "plasmid_labels.csv")
        if not os.path.exists(p):
            continue
        parts = t.split("_")
        prep = parts[1] if len(parts) > 2 else "?"
        depth = t.split("_d")[-1]
        iso = parts[0]
        for r in csv.DictReader(open(p)):
            g = truth.get((iso, r["plasmid"]))
            L = int(r["plasmid_len"])
            rows.append({
                "tag": t, "isolate": iso, "prep": prep, "depth": depth,
                "plasmid": r["plasmid"], "plasmid_len": L, "label": r["label"],
                "stratum": (g or {}).get("stratum", ""),
                "analysis_set": analysis_set_of(g, L),
                "cointegrated": (g or {}).get("cointegrated", ""),
                "repeat_role": (g or {}).get("repeat_role", ""),
                "in_cohort": ("_lig_d" in t or "_rap_d" in t),
            })
    return rows


def cell_failure_rates(rows):
    """failure % per (prep, depth) cell over interventional units."""
    per = collections.defaultdict(collections.Counter)
    for r in rows:
        if not r["in_cohort"] or r["analysis_set"] != INTERVENTIONAL:
            continue
        per[(r["prep"], r["depth"])][r["label"]] += 1
    out = {}
    for k, c in per.items():
        n = sum(c.values())
        fail = n - c.get("recovered", 0) - c.get("gray_50_95", 0)
        out[k] = {"n": n, "fail": fail, "rate": 100 * fail / n if n else 0.0,
                  "counts": dict(c)}
    return out


def discordance(rows):
    """Same definition as the 06_baselines.py adequacy gate, recomputed here so
    the report is standalone and the two numbers can be cross-checked."""
    idx = {}
    for r in rows:
        if r["in_cohort"] and r["analysis_set"] == INTERVENTIONAL:
            idx[(r["isolate"], r["plasmid"], r["prep"], r["depth"])] = r
    preps = sorted({k[2] for k in idx})
    depths = sorted({k[3] for k in idx}, key=lambda d: int(d))
    units = sorted({(k[0], k[1]) for k in idx})
    fail = lambda r: r["label"] in ("absent", "fragmented", "absorbed")
    cls, binr = set(), set()
    for iso, p in units:
        for d in depths:
            g = [idx[(iso, p, pr, d)] for pr in preps if (iso, p, pr, d) in idx]
            if len(g) >= 2:
                if len({x["label"] for x in g}) > 1:
                    cls.add((iso, p))
                if len({fail(x) for x in g}) > 1:
                    binr.add((iso, p))
        for pr in preps:
            g = [idx[(iso, p, pr, d)] for d in depths if (iso, p, pr, d) in idx]
            if len(g) >= 2:
                if len({x["label"] for x in g}) > 1:
                    cls.add((iso, p))
                if len({fail(x) for x in g}) > 1:
                    binr.add((iso, p))
    return {"units": len(units), "class": len(cls), "binary": len(binr)}


def main():
    truth = load_truth()
    rows = collect_labels(truth)

    # ------------------------------------------------------------------ 0
    section("0. REALISED CLASS COUNTS AGAINST FLOORS")
    inter = [r for r in rows if r["analysis_set"] == INTERVENTIONAL]
    all_c = collections.Counter(r["label"] for r in rows)
    int_c = collections.Counter(r["label"] for r in inter)
    print("  Floors apply to INTERVENTIONAL units: the <4 kb band drops in all")
    print("  four cells regardless of prep or depth, so counting it toward a")
    print("  mechanism floor would count constants as observations.")
    print()
    print("  %-12s %10s %14s %8s  %s"
          % ("class", "all units", "interventional", "floor", "verdict"))
    floor_status = {}
    for k in ("recovered", "fragmented", "absorbed", "absent", "gray_50_95"):
        fl = CLASS_FLOORS.get(k)
        n = int_c.get(k, 0)
        if k in UNGATED_CLASSES:
            verdict = "REPORTED, NOT GATED"
        elif fl is None:
            verdict = "-- no floor"
        elif n >= fl:
            verdict = "FLOOR MET"
        else:
            verdict = "FLOOR NOT MET -- NOT TESTABLE"
        if fl is not None:
            floor_status[k] = (n, fl, n >= fl)
        print("  %-12s %10d %14d %8s  %s"
              % (k, all_c.get(k, 0), n, fl if fl else "-", verdict))
    print()
    for k, why in UNGATED_CLASSES.items():
        print("  '%s' is REPORTED, NOT GATED: %s." % (k, why))
        print("     Its floor was withdrawn on probe evidence, not relaxed after a")
        print("     near miss. No per-class claim for '%s' rests on this arm." % k)
    print()
    for k, (n, fl, okc) in floor_status.items():
        if not okc:
            print("  '%s' has %d interventional units against a floor of %d."
                  % (k, n, fl))
            print("     NO per-class performance claim may be made for '%s'." % k)
    if floor_status and all(v[2] for v in floor_status.values()):
        print("  All remaining class floors met -- those per-class claims are admissible.")

    # ----------------------------------------------------------------- 0b
    section("0b. REBUILD GATE PANEL (all three must pass)")
    d = discordance(rows)
    cells = cell_failure_rates(rows)
    rates = [v["rate"] for v in cells.values()]
    spread = (max(rates) - min(rates)) if rates else 0.0
    gates = [
        ("discordant pairs >= %d" % MIN_DISCORDANT,
         d["class"] >= MIN_DISCORDANT,
         "%d discordant plasmids (binary flips %d) over %d interventional units"
         % (d["class"], d["binary"], d["units"])),
        ("absorbed units >= %d" % CLASS_FLOORS["absorbed"],
         int_c.get("absorbed", 0) >= CLASS_FLOORS["absorbed"],
         "%d interventional absorbed units" % int_c.get("absorbed", 0)),
        ("cell failure rates SPREAD (range >= %.0f pp)" % MIN_CELL_SPREAD_PP,
         spread >= MIN_CELL_SPREAD_PP,
         "range %.1f pp across %d cells (%s)"
         % (spread, len(rates), ", ".join("%.1f" % r for r in sorted(rates)))),
    ]
    for name, okg, detail in gates:
        print("  [%s] %-42s %s" % ("PASS" if okg else "FAIL", name, detail))
    n_pass = sum(1 for _, okg, _ in gates if okg)
    print()
    print("  'fragmented >= 20' was the fourth gate and has been WITHDRAWN on")
    print("  probe evidence (0/164 across three probes). It is not counted here.")
    print()
    if n_pass == len(gates):
        print("  ALL GATES PASS -- the cohort is adequate for what this arm claims:")
        print("  pipeline validation, and prep/depth effects on absence.")
        print("  It does NOT follow that B3/B4 are worth running. Check section 2:")
        print("  v2 held 45 real edges across 200 graphs with 175 edgeless, and a")
        print("  message-passing model on an edgeless corpus is an MLP. The")
        print("  relational-topology claim needs real assembly graphs.")
    else:
        print("  %d of %d GATES FAIL -- REBUILD THE COHORT." % (len(gates) - n_pass,
                                                               len(gates)))
        print("  Do not run B3/B4 on it: numbers from a cohort that fails these")
        print("  gates are discarded at the next rebuild anyway.")

    # ------------------------------------------------------------------ 1
    section("1. REALISED CLASS DISTRIBUTION (four label classes)")
    tot = sum(all_c.values())
    for k in ("recovered", "fragmented", "absorbed", "absent", "gray_50_95"):
        v = all_c.get(k, 0)
        print("  %-12s %4d  %5.1f%%" % (k, v, 100 * v / tot if tot else 0))
    print("  %-12s %4d" % ("TOTAL", tot))
    print()
    print("  by analysis set:")
    for aset in sorted({r["analysis_set"] for r in rows}):
        c = collections.Counter(r["label"] for r in rows
                                if r["analysis_set"] == aset)
        print("    %-20s n=%-5d %s" % (aset, sum(c.values()), dict(c)))
    print()
    print("  by stratum:")
    for s in ("tiny", "small", "medium", "large"):
        sub = [r for r in rows if r["stratum"] == s]
        if not sub:
            continue
        c = collections.Counter(r["label"] for r in sub)
        fail = sum(v for k, v in c.items() if k not in ("recovered", "gray_50_95"))
        print("    %-7s n=%-5d failure=%5.1f%%  %s"
              % (s, len(sub), 100 * fail / len(sub), dict(c)))
    print()
    print("  by cell (prep, depth) -- INTERVENTIONAL units only:")
    for k in sorted(cells):
        v = cells[k]
        print("    prep=%-4s depth=%-3sx  n=%-4d failure=%4.1f%%  %s"
              % (k[0], k[1], v["n"], v["rate"], v["counts"]))
    print("    range across cells: %.1f pp  (floor %.0f pp)"
          % (spread, MIN_CELL_SPREAD_PP))

    # ------------------------------------------------------------------ 1b
    section("1b. REPEAT-STRUCTURE YIELD (v2 induced mechanisms)")
    carriers = [r for r in rows if r["repeat_role"] == "frag"]
    coints = [r for r in rows if r["cointegrated"] == "true"]
    if not (carriers or coints):
        print("  no repeat structure in this cohort (ground_truth.csv carries no")
        print("  repeat_role / cointegrated columns -- this is a v1 cohort)")
    else:
        print("  IS-carrier units (role=frag)  n=%-4d %s"
              % (len(carriers),
                 dict(collections.Counter(r["label"] for r in carriers))))
        print("  cointegrated units            n=%-4d %s"
              % (len(coints),
                 dict(collections.Counter(r["label"] for r in coints))))
        if coints:
            hit = sum(1 for r in coints if r["label"] == "absorbed")
            print("  -> cointegrate -> absorbed conversion: %d/%d = %.0f%%"
                  % (hit, len(coints), 100 * hit / len(coints)))
        if carriers:
            hit = sum(1 for r in carriers if r["label"] == "fragmented")
            print("  -> IS carrier -> fragmented conversion: %d/%d = %.0f%%"
                  % (hit, len(carriers), 100 * hit / len(carriers)))

    # ------------------------------------------------------------------ 2
    section("2. NODE AND EDGE COUNTS PER GRAPH")
    diag = json.load(open(os.path.join(RES, "graph_diagnostics.json")))
    ns = [x["n_segments"] for x in diag]
    es = [x["n_simple_edges"] for x in diag]
    print("  graphs: %d" % len(diag))
    print("  segments  min=%d median=%d max=%d" % (min(ns), statistics.median(ns), max(ns)))
    print("  edges     min=%d median=%d max=%d" % (min(es), statistics.median(es), max(es)))
    print("  isolated nodes  median=%d max=%d"
          % (statistics.median([x["isolated_nodes"] for x in diag]),
             max(x["isolated_nodes"] for x in diag)))
    print()
    print("  %-20s %8s %7s %10s %6s %9s" %
          ("tag", "segments", "edges", "selfloops", "comps", "chromset"))
    for x in sorted(diag, key=lambda z: z["tag"]):
        print("  %-20s %8d %7d %10d %6d %9d%s"
              % (x["tag"], x["n_segments"], x["n_simple_edges"], x["self_loops"],
                 x["n_components"], x["chrom_set_size"],
                 "  FALLBACK" if x["chrom_fallback"] else ""))

    section("3. FRAGMENTED-CHROMOSOME FALLBACK")
    fb = [x["tag"] for x in diag if x["chrom_fallback"]]
    print("  graphs with no segment >=500 kb (fell back to all-contig median): %d"
          % len(fb))
    for t in fb:
        print("    ", t)
    if not fb:
        print("  none -- every graph had a chromosome-sized segment")

    section("4. SELF-LOOP COUNTS (diagnostic only, never a feature)")
    sl = [x["self_loops"] for x in diag]
    print("  total self-loops across corpus: %d" % sum(sl))
    print("  per graph  min=%d median=%d max=%d" % (min(sl), statistics.median(sl), max(sl)))
    print("  (a circular contig is represented by a self-loop; these are dropped")
    print("   from the degree feature so circularity cannot leak via the graph)")

    section("5. SIMULATED-VERSUS-REAL GRAPH STATISTICS")
    print("  PENDING -- requires real assemblies. No real cohort has been")
    print("  assembled yet (no ONT reads on disk for Wick / Johnson / NEKSUS).")
    print("  The comparison is defined and will run once the first real graphs exist.")
    print()
    print("  What IS available for a partial reference: the Wick deposited")
    print("  reference assemblies (7 isolates) give replicon-size distributions:")
    # NB the 10 kb cut below is Johnson's/Wick's reporting convention for
    # 'small plasmid', kept so this comparison stays commensurable with the
    # published tables. It is NOT the cohort's 'small' stratum bound, which is
    # 4-12 kb (see 02_build_sim_refs.py). Do not unify the two.
    tp = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\Dataset\assemblies.tar.gz"
    if os.path.exists(tp):
        t = tarfile.open(tp)
        sizes, small = [], 0
        for m in t.getmembers():
            fh = t.extractfile(m)
            ln, name = 0, None
            for raw in fh:
                line = raw.decode("utf-8", "replace").strip()
                if line.startswith(">"):
                    if name and ln < 500000:
                        sizes.append(ln)
                        if ln < 10000:
                            small += 1
                    name, ln = line, 0
                else:
                    ln += len(line)
            if name and ln < 500000:
                sizes.append(ln)
                if ln < 10000:
                    small += 1
        print("    Wick reference plasmids: n=%d, <10 kb: %d, median %d bp"
              % (len(sizes), small, statistics.median(sizes)))
    gt = list(truth.values())
    ss = [int(r["plasmid_len"]) for r in gt]
    if ss:
        print("    simulated plasmids     : n=%d, <10 kb: %d, median %d bp"
              % (len(ss), sum(1 for x in ss if x < 10000), statistics.median(ss)))
        band = collections.Counter(r.get("stratum", "?") for r in gt)
        print("    simulated by stratum   :", dict(band))

    section("6. PER-MODE FLYE RECOVERY (correctness gate)")
    print("  PENDING -- the gate is defined against Johnson's 14 isolates at full")
    print("  depth, whose ONT reads are not on disk. Gate target, now pinned:")
    print("     flye-raw, replicate 1, 30 of 33 plasmids, tolerance +/-1")
    print("     (Johnson across-replicate range for this mode: 29-31)")
    print()
    print("  Both Flye modes ARE already exercised in simulation:")
    mo = collections.Counter()
    for dd in sorted(os.listdir(os.path.join(ROOT, "asm"))):
        mp = os.path.join(ROOT, "asm", dd, "manifest.json")
        if os.path.exists(mp):
            m = json.load(open(mp))
            mo[m.get("flye_mode", "?")] += 1
    for k, v in sorted(mo.items()):
        print("     %-12s %d assemblies" % (k, v))

    section("7. LABEL VALIDATION AGAINST GROUND TRUTH")
    lv = os.path.join(RES, "label_validation.csv")
    if os.path.exists(lv):
        lrows = list(csv.DictReader(open(lv)))
        by = collections.defaultdict(list)
        for r in lrows:
            by[r["label"]].append(float(r["realised_rel_depth"]))
        print("  realised relative depth by assigned label:")
        for lab in sorted(by):
            v = sorted(by[lab])
            print("    %-12s n=%-4d median=%8.2f" % (lab, len(v), v[len(v) // 2]))
    else:
        print("  not run yet")

    section("8. BASELINE LADDER")
    bp = os.path.join(RES, "baselines.json")
    if os.path.exists(bp):
        b = json.load(open(bp))
        if b.get("gate") == "FAILED":
            print("  NOT RUN -- the adequacy gate in 06_baselines.py refused this")
            print("  cohort: %d discordant plasmids against a floor of %d."
                  % (b["discordance"]["discordant_class"],
                     b.get("min_discordant", MIN_DISCORDANT)))
        else:
            dd = b.get("dataset", {})
            print("  n = %s ISOLATES (resampling unit)" % dd.get("n_isolates", "?"))
            print("  nested interventional units=%s  positives=%s  grey held out=%s"
                  % (dd.get("n_units_nested_interventional"),
                     dd.get("n_positive_interventional"),
                     dd.get("gray_band_held_out")))
            for arm in ("with_is_empty", "without_is_empty"):
                if not b.get(arm):
                    continue
                print("  %s (PR-AUC, leave-isolate-out):" % arm)
                for k, v in b[arm]["pr_auc"].items():
                    print("     %-16s %.3f" % (k, v))
            stab = b.get("cross_arm_stability") or {}
            if stab:
                print("  cross-arm stability:")
                for k, v in stab.items():
                    print("     %-12s %s" % (k, v["verdict"]))
                unres = [k for k, v in stab.items() if v["verdict"] == "UNRESOLVED"]
                if unres:
                    print("     UNRESOLVED comparisons must not be reported as")
                    print("     clearing or failing: %s" % ", ".join(unres))
    else:
        print("  not run yet")


if __name__ == "__main__":
    main()
