"""
v3 DESIGN PROBE, ROUND 2 -- why did nothing fragment?

Round 1 result: two copies of a repeat at 3, 8, 15 and 25 kb, at 99 % identity,
in 94-170 kb carriers, ALL came back `recovered` with a single support segment.
Flye assembled every carrier complete and circular. Repeat LENGTH alone does
not induce fragmentation.

Round 1 tested the case most favourable to the assembler and least favourable
to the mechanism, on four axes at once. This round moves each one, so a second
negative means something rather than nothing:

  IDENTITY. 99 % over a 15 kb repeat is 150 differences -- ample signal for
  Flye to separate the copies, and more divergence than the reads carry
  (nanopore2023 is ~1 % error). Recent IS transposition gives copies that are
  100 % identical, and an exactly identical repeat cannot be separated by
  sequence at all, only by coverage. This is the axis I expect to matter most.

  COPY COUNT. Two copies in a circle is the easiest possible motif: the repeat
  sits at 2x coverage, the two unique arms at 1x, and the ratio resolves it.
  Three or more copies makes the coverage argument ambiguous.

  CROSS-REPLICON PLACEMENT. Round 1 lifted the family from the chromosome, so a
  native locus existed, but added no extra chromosomal copies. v2's design put
  2-4 in the chromosome. Ambiguity between a plasmid copy and a chromosome copy
  is the thing that actually strands a plasmid contig.

  CHEMISTRY AND DEPTH. Round 1 ran r10_hq at 30x -- the cleanest reads and the
  most coverage, i.e. the easiest resolution. The cohort is half r9_raw, and
  15x is a rung. probe_v3b.sh runs r9_raw at 15x as the second cell.

If nothing fragments here either, the frag mechanism is not available at this
read length and the honest move is to re-scope, not to keep tuning.

Writes to work/probe2/ -- both the v2 corpus and round 1 are left untouched.
"""
import os, sys, csv, random, json

ROOT      = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\work"
POOL      = os.path.join(ROOT, "refs", "plasmid_pool.fasta")
ECOLI_DIR = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\Dataset\ecoli_10\ncbi_dataset\data"
PROBE     = os.path.join(ROOT, "probe2")
SIM       = os.path.join(PROBE, "sim")

SEED = 20260815
BASES = "ACGT"

# name, repeat_len, copies_in_plasmid, copies_in_chromosome, identity
# Only ONE axis moves between neighbouring rungs wherever possible.
LADDER_C = [
    ("c1_len15_n2_id990",  15_000, 2, 2, 0.990),   # round-1 setting + chrom copies
    ("c2_len15_n2_id1000", 15_000, 2, 2, 1.000),   # <- identity is the only change
    ("c3_len15_n3_id1000", 15_000, 3, 2, 1.000),   # <- copy count is the only change
    ("c4_len25_n2_id1000", 25_000, 2, 2, 1.000),   # <- length is the only change
    ("c5_len05_n4_id1000",  5_000, 4, 3, 1.000),   # short but many, IS-like
    ("c6_len02_n6_id1000",  2_000, 6, 4, 1.000),   # true IS length, many copies
]
CARRIER_MIN_MULT = 3      # carrier must be >= this x (repeat_len * copies)


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def mutate(seq, rng, identity):
    """Substitutions only. identity=1.0 returns the sequence untouched, which
    is the whole point of this round -- an exactly identical repeat."""
    if identity >= 1.0:
        return seq
    s = list(seq)
    n = int(round(len(s) * (1.0 - identity)))
    for p in rng.sample(range(len(s)), min(n, len(s))):
        cur = s[p].upper()
        s[p] = rng.choice([b for b in BASES if b != cur]) if cur in BASES else cur
    return "".join(s)


def insert_all(seq, items):
    """Insert every copy in ONE descending pass so earlier offsets stay valid."""
    out, rec = seq, []
    for p, cp in sorted(items, key=lambda x: -x[0]):
        out = out[:p] + cp + out[p:]
        rec.append({"at": p, "len": len(cp)})
    return out, rec


def spread(rng, seq_len, n, margin=3000):
    """Evenly spaced, jittered, off the ends -- every copy flanked by unique
    sequence so it is a genuine dispersed repeat, not a tandem array."""
    span = seq_len - 2 * margin
    step = span // (n + 1)
    return [margin + step * (k + 1) + rng.randint(-step // 4, step // 4)
            for k in range(n)]


def main():
    rng = random.Random(SEED)
    os.makedirs(SIM, exist_ok=True)

    chrom = None
    for d in sorted(os.listdir(ECOLI_DIR)):
        sub = os.path.join(ECOLI_DIR, d)
        if not os.path.isdir(sub):
            continue
        for f in sorted(os.listdir(sub)):
            if f.endswith((".fna", ".fasta")):
                recs = sorted(read_fasta(os.path.join(sub, f)), key=lambda x: -len(x[1]))
                if recs and len(recs[0][1]) > 4_000_000:
                    chrom = (d, recs[0][1])
                break
        if chrom:
            break
    if not chrom:
        sys.exit("no chromosome backbone found")
    cacc, cseq = chrom
    print("chromosome: %s  %d bp" % (cacc, len(cseq)))

    pool = sorted(((h.split()[0], s) for h, s in read_fasta(POOL)),
                  key=lambda x: -len(x[1]))
    used = set()
    truth, placements = [], []
    plasmid_records = []
    chrom_items = []          # ALL chromosome insertions, applied in one pass

    print()
    print("LADDER C -- identity / copy count / cross-replicon placement")
    for name, rlen, npl, nch, ident in LADDER_C:
        need = rlen * npl * CARRIER_MIN_MULT
        cand = [(a, s) for a, s in pool if a not in used and len(s) >= need]
        if not cand:
            print("  %-20s NO CARRIER >= %d bp -- skipped" % (name, need))
            continue
        acc, seq = cand[-1]                    # smallest carrier that still fits
        used.add(acc)

        src = rng.randrange(0, len(cseq) - rlen)
        fam = cseq[src:src + rlen]

        pos = spread(rng, len(seq), npl)
        newseq, rec = insert_all(seq, [(p, mutate(fam, rng, ident)) for p in pos])
        plasmid_records.append((name, newseq))

        cpos = spread(rng, len(cseq), nch, margin=50_000)
        for p in cpos:
            chrom_items.append((p, mutate(fam, rng, ident)))

        truth.append({"plasmid_name": name, "plasmid_acc": acc,
                      "plasmid_len": len(newseq), "carrier_len": len(seq),
                      "repeat_len": rlen, "copies_in_plasmid": npl,
                      "copies_in_chrom": nch, "identity": ident,
                      "copy_number": 2})
        placements.append({"plasmid": name, "family_src": src, "repeat_len": rlen,
                           "identity": ident, "in_plasmid": rec,
                           "in_chrom": sorted(cpos)})
        print("  %-20s repeat %5d bp x%d in plasmid, x%d in chrom, id=%.3f  "
              "carrier %d -> %d bp" % (name, rlen, npl, nch, ident, len(seq), len(newseq)))

    # chromosome: ONE insertion pass for every family at once
    newc, crec = insert_all(cseq, chrom_items)
    print()
    print("chromosome: %d insertions, %d -> %d bp" % (len(crec), len(cseq), len(newc)))

    ref = os.path.join(SIM, "probe02_ref.fasta")
    with open(ref, "w") as fo:
        fo.write(">chromosome depth=1.0 circular=true\n")
        for k in range(0, len(newc), 70):
            fo.write(newc[k:k + 70] + "\n")
        for name, seq in plasmid_records:
            fo.write(">%s depth=2 circular=true\n" % name)
            for k in range(0, len(seq), 70):
                fo.write(seq[k:k + 70] + "\n")

    back = {h.split()[0]: s for h, s in read_fasta(ref)}
    bad = [t["plasmid_name"] for t in truth
           if len(back.get(t["plasmid_name"], "")) != t["plasmid_len"]]
    if bad or len(back) != len(truth) + 1:
        sys.exit("VERIFY FAILED: %r  (records=%d)" % (bad, len(back)))
    # an identity=1.000 copy must be findable verbatim in the written plasmid
    for p in placements:
        if p["identity"] < 1.0:
            continue
        fam = cseq[p["family_src"]:p["family_src"] + p["repeat_len"]]
        if back[p["plasmid"]].count(fam) < 1:
            sys.exit("VERIFY FAILED: exact repeat not found in %s" % p["plasmid"])
    print("verified: identical copies present verbatim where identity=1.000")

    with open(os.path.join(SIM, "probe_truth.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(truth[0].keys()))
        w.writeheader()
        w.writerows(truth)
    json.dump({"seed": SEED, "chromosome": cacc, "chromosome_len": len(newc),
               "placements": placements},
              open(os.path.join(SIM, "probe_design.json"), "w"), indent=1)

    total = len(newc) + sum(t["plasmid_len"] for t in truth)
    print("records %d, %.2f Mb reference" % (len(back), total / 1e6))
    print("-> %s" % ref)


if __name__ == "__main__":
    main()
