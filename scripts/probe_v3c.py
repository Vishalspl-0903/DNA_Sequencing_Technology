"""
v3 DESIGN PROBE, ROUND 3 -- map the prep window in TWO dimensions.

Round 1 (ladder A) pinned copy number at 15x and found the label flipping on
prep between 9,083 and 10,046 bp. That is the right measurement for isolating
length, but it is not enough to set a stratum, because the cohort draws copy
number at random over 8-30x and DEPTH is what the ligation bias actually acts
on. A 9 kb plasmid at 30x and a 9 kb plasmid at 8x are not the same experiment.

That interaction is precisely what diluted cohort v2: the 'small' stratum was
4-12 kb x copy 8-30, so it mixed plasmids that fail in both preps, plasmids
that pass in both, and the thin band that actually flips -- giving a 5.5 pp
cell spread against a 10 pp floor.

This maps length x copy-number over the FULL 2x2 the gate is computed on
(lig/rap x 15x/30x), so STRATUM_BOUNDS and COPY_RANGE can be set from measured
values rather than from one slice.

Plasmids span ~7.5-12.5 kb with copy number assigned round-robin over
{8, 15, 25}, so length and copy number are crossed rather than confounded.

Writes to work/probe3/ -- v2, round 1 and round 2 are all untouched.
"""
import os, sys, csv, json

ROOT      = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\work"
POOL      = os.path.join(ROOT, "refs", "plasmid_pool.fasta")
ECOLI_DIR = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\Dataset\ecoli_10\ncbi_dataset\data"
PROBE     = os.path.join(ROOT, "probe3")
SIM       = os.path.join(PROBE, "sim")

LEN_LO, LEN_HI = 7_500, 12_500
COPIES = [8, 15, 25]          # spans the cohort's COPY_RANGE for 'small'
N_MAX = 27                    # 9 per copy level


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


def main():
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

    cand = sorted(((a.split()[0], s) for a, s in read_fasta(POOL)
                   if LEN_LO <= len(s) <= LEN_HI), key=lambda x: len(x[1]))
    if len(cand) < len(COPIES) * 3:
        sys.exit("only %d pool plasmids in %d-%d bp" % (len(cand), LEN_LO, LEN_HI))

    # thin evenly across the length range, then assign copy number round-robin
    # so length and copy number are CROSSED, not confounded
    step = max(1, len(cand) // N_MAX)
    picked = cand[::step][:N_MAX]

    truth, records = [], [(">chromosome depth=1.0 circular=true", cseq)]
    print("chromosome: %s  %d bp" % (cacc, len(cseq)))
    print("pool candidates in %d-%d bp: %d, taking %d"
          % (LEN_LO, LEN_HI, len(cand), len(picked)))
    print()
    print("  %-22s %8s %6s" % ("name", "len", "copy"))
    for i, (acc, seq) in enumerate(picked):
        cp = COPIES[i % len(COPIES)]
        name = "p%02d_len%05d_cp%02d" % (i, len(seq), cp)
        records.append((">%s depth=%d circular=true" % (name, cp), seq))
        truth.append({"plasmid_name": name, "plasmid_acc": acc,
                      "plasmid_len": len(seq), "copy_number": cp})
        print("  %-22s %8d %6d" % (name, len(seq), cp))

    ref = os.path.join(SIM, "probe03_ref.fasta")
    with open(ref, "w") as fo:
        for hdr, seq in records:
            fo.write(hdr + "\n")
            for k in range(0, len(seq), 70):
                fo.write(seq[k:k + 70] + "\n")

    back = {h.split()[0]: s for h, s in read_fasta(ref)}
    bad = [t["plasmid_name"] for t in truth
           if len(back.get(t["plasmid_name"], "")) != t["plasmid_len"]]
    if bad or len(back) != len(records):
        sys.exit("VERIFY FAILED: %r" % (bad,))

    with open(os.path.join(SIM, "probe_truth.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(truth[0].keys()))
        w.writeheader()
        w.writerows(truth)
    json.dump({"chromosome": cacc, "chromosome_len": len(cseq),
               "len_range": [LEN_LO, LEN_HI], "copies": COPIES},
              open(os.path.join(SIM, "probe_design.json"), "w"), indent=1)

    total = len(cseq) + sum(t["plasmid_len"] for t in truth)
    print()
    print("verified: %d records, %d plasmids, %.2f Mb" % (len(back), len(truth), total / 1e6))
    print("copy levels:", {c: sum(1 for t in truth if t["copy_number"] == c) for c in COPIES})
    print("-> %s" % ref)


if __name__ == "__main__":
    main()
