"""
Unit test for the labeller (step 4), run before any real graph is trusted.

Builds four synthetic cases with known answers and checks the label that comes
back. This validates the parts most likely to be silently wrong: merged-interval
coverage, the absorbed-vs-recovered precedence, and the support-set thresholds.

  case A  recovered  : plasmid sits on one segment of its own size
  case B  fragmented : plasmid split across two segments, ~50/50
  case C  absorbed   : plasmid embedded in a segment >=3x its length
  case D  absent     : plasmid not in the graph at all
"""
import os, sys, random, shutil, importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("lab", os.path.join(HERE, "04_label.py"))
lab = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lab)

TMP = os.environ.get("TEST_TMP", "/tmp/labeller_test")
MM2 = sys.argv[1] if len(sys.argv) > 1 else "minimap2"

rng = random.Random(7)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(s, rate=0.01):
    out = []
    for c in s:
        if rng.random() < rate:
            out.append(rng.choice("ACGT"))
        else:
            out.append(c)
    return "".join(out)


def wrap(s, w=70):
    return "\n".join(s[i:i + w] for i in range(0, len(s), w))


def main():
    if os.path.exists(TMP):
        shutil.rmtree(TMP)
    asm = os.path.join(TMP, "asm")
    os.makedirs(asm)

    pA = seq(6000)                      # recovered
    pB = seq(8000)                      # fragmented
    pC = seq(4000)                      # absorbed
    pD = seq(5000)                      # absent

    filler = seq(30000)
    segments = {
        # own-size segment -> recovered
        "edge_1": mutate(pA),
        # two halves -> fragmented (overlap-free split)
        "edge_2": mutate(pB[:4200]),
        "edge_3": mutate(pB[3900:]),
        # 4 kb plasmid inside a 30 kb segment -> ratio 7.5x -> absorbed
        "edge_4": mutate(filler[:13000] + pC + filler[13000:]),
        # unrelated
        "edge_5": seq(20000),
    }

    with open(os.path.join(asm, "assembly_graph.gfa"), "w") as fo:
        for n, s in segments.items():
            fo.write("S\t%s\t%s\tdp:i:40\n" % (n, s))
        fo.write("L\tedge_2\t+\tedge_3\t+\t0M\n")
        fo.write("L\tedge_1\t+\tedge_1\t+\t0M\n")     # self-loop, must be dropped

    ref = os.path.join(TMP, "sim_ref.fasta")
    with open(ref, "w") as fo:
        fo.write(">chromosome depth=1.0 circular=true\n%s\n" % wrap(seq(600000)))
        for name, s in (("plasmid_1", pA), ("plasmid_2", pB),
                        ("plasmid_3", pC), ("plasmid_4", pD)):
            fo.write(">%s depth=10 circular=true\n%s\n" % (name, wrap(s)))

    out = os.path.join(TMP, "out")
    summary = lab.label_one(asm, ref, out, MM2)
    if summary is None:
        print("FAIL: labeller returned nothing"); return 1

    import csv
    got = {r["plasmid"]: r for r in
           csv.DictReader(open(os.path.join(out, "plasmid_labels.csv")))}

    expect = {"plasmid_1": "recovered", "plasmid_2": "fragmented",
              "plasmid_3": "absorbed", "plasmid_4": "absent"}

    print("=" * 64)
    print("LABELLER UNIT TEST")
    print("=" * 64)
    ok = True
    for p in sorted(expect):
        g = got.get(p, {})
        lab_got = g.get("label", "<missing>")
        good = (lab_got == expect[p])
        ok &= good
        print("  %-10s expect=%-11s got=%-11s best=%.2f cum=%.2f nsup=%s  %s"
              % (p, expect[p], lab_got,
                 float(g.get("best_cov_frac", 0)),
                 float(g.get("cumulative_cov_frac", 0)),
                 g.get("n_support", "?"),
                 "OK" if good else "*** MISMATCH ***"))
    print()
    print("  graph-head count (absent):", summary["n_missing_for_graph_head"],
          "expect 1")
    ok &= (summary["n_missing_for_graph_head"] == 1)

    # node labels: only fragmented/absorbed segments are positive
    nodes = {r["segment"]: r for r in
             csv.DictReader(open(os.path.join(out, "node_labels.csv")))}
    pos = sorted(n for n, r in nodes.items() if r["node_positive"] == "1")
    print("  positive nodes:", pos, "expect ['edge_2','edge_3','edge_4']")
    ok &= (pos == ["edge_2", "edge_3", "edge_4"])

    print()
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
