"""
v3 DESIGN PROBE -- read out the two ladders.

Answers exactly two questions, and refuses to answer anything else:

  A. At FRAGLEN 15000,13000, over what plasmid length range do ligation and
     rapid actually produce DIFFERENT labels? That range, and only that range,
     is where the interventional stratum belongs. v2 put 4-12 kb there and got
     a 5.5 pp cell spread against a 10 pp floor.

  B. How long must a duplicated repeat be before Flye stops resolving it and
     the carrier plasmid fragments? v2 used 1.0-1.9 kb elements against 15 kb
     reads and converted 0/148. This also decides whether the corpus has any
     edges at all: v2 had 45 across 200 graphs, 91 % of nodes isolated.

Both ladders vary ONE factor. Ladder A holds copy number fixed at 15x so length
is the only variable; ladder B holds carrier band, copy number and identity
fixed so repeat length is the only variable.
"""
import os, csv, json, collections

ROOT  = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\work"
PROBE = os.path.join(ROOT, "probe")
CELLS = ["probe01_lig_d30", "probe01_rap_d30"]
PREP  = {"probe01_lig_d30": "lig", "probe01_rap_d30": "rap"}


def load_labels(tag):
    p = os.path.join(PROBE, "graphs", tag, "plasmid_labels.csv")
    if not os.path.exists(p):
        return None
    return {r["plasmid"]: r for r in csv.DictReader(open(p))}


def load_yield(tag):
    """Realised bases per replicon, straight from the read source tags.

    This is the prep effect at the DEPTH level, before the assembler gets a
    say. v1 measured 17x depletion here that never reached the label level --
    which is exactly the failure mode this probe is checking for.
    """
    p = os.path.join(PROBE, "asm", tag, "replicon_yield.tsv")
    if not os.path.exists(p):
        return {}
    out = {}
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3:
            out[f[0]] = int(f[2])
    return out


def gfa_stats(tag):
    p = os.path.join(PROBE, "asm", tag, "assembly_graph.gfa")
    if not os.path.exists(p):
        return None
    seg = links = selfloops = 0
    for line in open(p):
        if line.startswith("S\t"):
            seg += 1
        elif line.startswith("L\t"):
            f = line.split("\t")
            if len(f) > 3 and f[1] == f[3]:
                selfloops += 1
            else:
                links += 1
    return {"segments": seg, "links": links, "self_loops": selfloops}


def main():
    truth = list(csv.DictReader(open(os.path.join(PROBE, "sim", "probe_truth.csv"))))
    design = json.load(open(os.path.join(PROBE, "sim", "probe_design.json")))
    clen = design["chromosome_len"]

    labels = {t: load_labels(t) for t in CELLS}
    yields = {t: load_yield(t) for t in CELLS}
    gfas = {t: gfa_stats(t) for t in CELLS}

    missing = [t for t in CELLS if labels[t] is None]
    if missing:
        raise SystemExit("no labels for %r -- run probe_v3.sh then probe_v3_label.sh"
                         % missing)

    print("=" * 78)
    print("GRAPH SHAPE")
    print("=" * 78)
    print("  %-22s %9s %7s %11s" % ("cell", "segments", "links", "self-loops"))
    for t in CELLS:
        g = gfas[t]
        print("  %-22s %9d %7d %11d" % (t, g["segments"], g["links"], g["self_loops"]))
    print()
    print("  v2 cohort for comparison: 45 links across 200 graphs, 175 of them")
    print("  with zero. Links are what a message-passing model has to work with.")

    def rel_depth(t, name, plen):
        y = yields[t]
        cb = y.get("chromosome", 0)
        if not cb or not plen:
            return float("nan")
        return (y.get(name, 0) / plen) / (cb / clen)

    # ------------------------------------------------------------ ladder A
    print()
    print("=" * 78)
    print("LADDER A -- WHERE DOES PREP CHANGE THE LABEL?   (copy number fixed 15x)")
    print("=" * 78)
    print("  %-7s %8s | %8s %8s | %-11s %-11s %s"
          % ("rung", "len", "dep lig", "dep rap", "lig", "rap", ""))
    a = [t for t in truth if t["ladder"] == "A"]
    a.sort(key=lambda r: int(r["plasmid_len"]))
    split_lo = split_hi = None
    for r in a:
        n, L = r["plasmid_name"], int(r["plasmid_len"])
        ll = labels[CELLS[0]].get(n, {}).get("label", "?")
        lr = labels[CELLS[1]].get(n, {}).get("label", "?")
        fl = ll in ("absent", "fragmented", "absorbed")
        fr = lr in ("absent", "fragmented", "absorbed")
        mark = ""
        if fl != fr:
            mark = "<== PREP SPLIT"
            split_lo = L if split_lo is None else min(split_lo, L)
            split_hi = L if split_hi is None else max(split_hi, L)
        elif fl and fr:
            mark = "both fail"
        else:
            mark = "both ok"
        print("  %-7s %8d | %8.1f %8.1f | %-11s %-11s %s"
              % (n.replace("ladderA_", ""), L,
                 rel_depth(CELLS[0], n, L), rel_depth(CELLS[1], n, L), ll, lr, mark))

    print()
    if split_lo is None:
        print("  NO rung splits on prep. The interventional stratum cannot be")
        print("  placed from this probe -- widen ladder A or change FRAGLEN.")
    else:
        print("  PREP SPLIT SPANS %d - %d bp." % (split_lo, split_hi))
        print("  That is where the interventional stratum belongs. Rungs that")
        print("  fail in BOTH preps are constants; so are rungs that pass in both.")

    # ------------------------------------------------------------ ladder B
    print()
    print("=" * 78)
    print("LADDER B -- HOW LONG MUST A REPEAT BE TO FRAGMENT?  (2 copies, 99%% id)")
    print("=" * 78)
    print("  %-10s %9s | %-12s %-12s %s"
          % ("repeat", "carrier", "lig", "rap", "support segments (lig/rap)"))
    b = [t for t in truth if t["ladder"] == "B"]
    b.sort(key=lambda r: int(r["repeat_len"]))
    first_frag = None
    for r in b:
        n = r["plasmid_name"]
        L = int(r["plasmid_len"])
        rl = int(r["repeat_len"])
        rowl = labels[CELLS[0]].get(n, {})
        rowr = labels[CELLS[1]].get(n, {})
        ll, lr = rowl.get("label", "?"), rowr.get("label", "?")
        sl, sr = rowl.get("n_support", "?"), rowr.get("n_support", "?")
        if (ll == "fragmented" or lr == "fragmented") and first_frag is None:
            first_frag = rl
        print("  %-10s %9d | %-12s %-12s %s / %s%s"
              % ("%d kb" % (rl // 1000), L, ll, lr, sl, sr,
                 "   <== FRAGMENTS" if "fragmented" in (ll, lr) else ""))

    print()
    if first_frag is None:
        print("  NOTHING FRAGMENTED, up to a %d kb repeat."
              % (max(int(r["repeat_len"]) for r in b) // 1000))
        print("  Repeat length is not sufficient to induce fragmentation at this")
        print("  read length. The frag mechanism needs a different lever --")
        print("  lower FRAGLEN, or drop the mechanism and re-scope the claim.")
    else:
        print("  FRAGMENTATION STARTS AT %d kb repeats." % (first_frag // 1000))
        print("  The frag role must use elements at or above this length, in")
        print("  carriers long enough to hold 2 copies and stay in band.")

    print()
    print("=" * 78)
    print("VERDICT FOR 02_build_sim_refs.py")
    print("=" * 78)
    if split_lo is not None:
        print("  interventional stratum -> %d - %d bp" % (split_lo, split_hi))
    else:
        print("  interventional stratum -> UNDETERMINED")
    if first_frag is not None:
        print("  frag element length    -> >= %d bp" % first_frag)
        print("  frag carrier           -> >= %d bp to hold 2 copies in band"
              % (3 * first_frag))
    else:
        print("  frag element length    -> NOT ACHIEVABLE at FRAGLEN 15000,13000")


if __name__ == "__main__":
    main()
