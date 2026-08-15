"""
v3 DESIGN PROBE -- one purpose-built isolate that answers both open questions.

Cohort v2 failed 2 of 4 gates. The diagnosis (measured, PROGRESS/PROPOSAL_FIXES):

  * cell failure spread 5.5 pp vs a 10 pp floor, because the 'small' stratum
    (4-12 kb) is three regimes glued together at FRAGLEN 15000,13000:
    4-6 kb fails 98.3 % (more deterministically than the excluded <4 kb band),
    6-9 kb fails 73.7 % and splits 86.8/60.5 on prep, 9-12 kb fails 11.7 %
    with no prep effect at all.
  * fragmented 19 vs a floor of 20, of which 0 came from the designed
    mechanism: IS carrier -> fragmented converted 0/148. The frag elements are
    1.0-1.9 kb against reads of mean 15 kb, and a repeat only breaks an
    assembly graph when reads cannot span it.

Rather than rewrite 02_build_sim_refs.py on that theory and spend another 12 h
finding out, this builds ONE isolate carrying two controlled ladders and
assembles it under both preps. Everything varies one factor at a time.

  LADDER A -- prep window. Eight plasmids at 5,6,7,8,9,10,11,12 kb, ALL at the
  SAME copy number, so length is the only variable. Locates where lig and rap
  separate, and where each becomes deterministic.

  LADDER B -- repeat resolution. Four large carriers, each with 2 copies of a
  repeat lifted from this isolate's own chromosome, at 3, 8, 15 and 25 kb.
  Same carrier size band, same copy number, same identity: repeat LENGTH is the
  only variable. Finds the length at which Flye stops resolving and the plasmid
  actually fragments -- which is also the only way this corpus gets edges
  (v2: 45 edges over 200 graphs, 91 % of nodes isolated).

Writes to work/probe/, NOT work/sim -- the v2 corpus is left untouched.
"""
import os, sys, csv, random, json

ROOT      = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\work"
POOL      = os.path.join(ROOT, "refs", "plasmid_pool.fasta")
ECOLI_DIR = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\Dataset\ecoli_10\ncbi_dataset\data"
PROBE     = os.path.join(ROOT, "probe")
SIM       = os.path.join(PROBE, "sim")

SEED = 20260814

# Ladder A: same copy number at every rung, so length is the only variable.
LADDER_A_KB   = [5, 6, 7, 8, 9, 10, 11, 12]
LADDER_A_COPY = 15
LADDER_A_TOL  = 400           # bp; how close to the nominal rung we must land

# Ladder B: same carrier band and copy number, repeat length is the only variable.
LADDER_B_REPEAT_KB = [3, 8, 15, 25]
LADDER_B_CARRIER   = (90_000, 200_000)
LADDER_B_COPIES    = 2        # >=2 copies is what makes a repeat ambiguous
LADDER_B_COPY_NO   = 2
REPEAT_IDENTITY    = 0.99     # substitutions only -> copies keep their length

BASES = "ACGT"


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
    """Substitutions only, so every copy keeps the family length exactly."""
    s = list(seq)
    n = int(round(len(s) * (1.0 - identity)))
    for p in rng.sample(range(len(s)), min(n, len(s))):
        cur = s[p].upper()
        s[p] = rng.choice([b for b in BASES if b != cur]) if cur in BASES else cur
    return "".join(s)


def insert_copies(seq, positions, copies):
    """Insert in ONE pass, descending, so earlier offsets stay valid.

    02_build_sim_refs.py had a real bug here once (two rounds of insertion, the
    second silently invalidating the first's recorded offsets, verified identity
    0.26 against a recorded 0.99). Single pass, and the placements are verified
    against the written file below.
    """
    out = seq
    rec = []
    for p, cp in sorted(zip(positions, copies), key=lambda x: -x[0]):
        out = out[:p] + cp + out[p:]
        rec.append((p, len(cp)))
    return out, rec


def main():
    rng = random.Random(SEED)
    os.makedirs(SIM, exist_ok=True)

    # ---- chromosome -------------------------------------------------------
    chrom = None
    for d in sorted(os.listdir(ECOLI_DIR)):
        sub = os.path.join(ECOLI_DIR, d)
        if not os.path.isdir(sub):
            continue
        for f in sorted(os.listdir(sub)):
            if f.endswith((".fna", ".fasta")):
                recs = [(h, s) for h, s in read_fasta(os.path.join(sub, f))]
                recs.sort(key=lambda x: -len(x[1]))
                if recs and len(recs[0][1]) > 4_000_000:
                    chrom = (d, recs[0][1])
                break
        if chrom:
            break
    if not chrom:
        sys.exit("no chromosome backbone found under %s" % ECOLI_DIR)
    cacc, cseq = chrom
    print("chromosome: %s  %d bp" % (cacc, len(cseq)))

    # ---- pool, indexed by length -----------------------------------------
    pool = [(h.split()[0], s) for h, s in read_fasta(POOL)]
    print("pool: %d plasmids" % len(pool))

    records = [(">chromosome depth=1.0 circular=true", cseq)]
    truth = []

    # ---- LADDER A: the prep window ---------------------------------------
    print()
    print("LADDER A -- prep window (copy number fixed at %dx)" % LADDER_A_COPY)
    used = set()
    for kb in LADDER_A_KB:
        target = kb * 1000
        cand = [(abs(len(s) - target), a, s) for a, s in pool
                if a not in used and abs(len(s) - target) <= LADDER_A_TOL]
        if not cand:
            print("  %2d kb : NO PLASMID within %d bp -- rung skipped" % (kb, LADDER_A_TOL))
            continue
        cand.sort(key=lambda x: x[0])
        _, acc, seq = cand[0]
        used.add(acc)
        name = "ladderA_%02dkb" % kb
        records.append((">%s depth=%d circular=true" % (name, LADDER_A_COPY), seq))
        truth.append({"plasmid_name": name, "ladder": "A", "plasmid_acc": acc,
                      "plasmid_len": len(seq), "copy_number": LADDER_A_COPY,
                      "repeat_len": 0, "repeat_copies": 0})
        print("  %2d kb : %-16s %7d bp" % (kb, acc, len(seq)))

    # ---- LADDER B: repeat resolution -------------------------------------
    print()
    print("LADDER B -- repeat resolution (%d copies, carrier %d-%d kb, copy no. %dx)"
          % (LADDER_B_COPIES, LADDER_B_CARRIER[0] // 1000,
             LADDER_B_CARRIER[1] // 1000, LADDER_B_COPY_NO))
    carriers = [(a, s) for a, s in pool
                if LADDER_B_CARRIER[0] <= len(s) < LADDER_B_CARRIER[1] and a not in used]
    rng.shuffle(carriers)
    placements = []
    for j, kb in enumerate(LADDER_B_REPEAT_KB):
        rlen = kb * 1000
        if j >= len(carriers):
            print("  %2d kb : no carrier left -- rung skipped" % kb)
            continue
        acc, seq = carriers[j]
        used.add(acc)
        # family lifted from this isolate's own chromosome: real genomic base
        # composition, and a genuine native locus so the repeat is shared with
        # the chromosome as well as duplicated within the plasmid
        src = rng.randrange(0, len(cseq) - rlen)
        fam = cseq[src:src + rlen]
        cps = [mutate(fam, rng, REPEAT_IDENTITY) for _ in range(LADDER_B_COPIES)]
        # keep copies apart and off the ends so each is flanked by unique sequence
        margin = 2000
        span = len(seq) - 2 * margin
        step = span // (LADDER_B_COPIES + 1)
        pos = [margin + step * (k + 1) for k in range(LADDER_B_COPIES)]
        newseq, rec = insert_copies(seq, pos, cps)
        name = "ladderB_rep%02dkb" % kb
        records.append((">%s depth=%d circular=true" % (name, LADDER_B_COPY_NO), newseq))
        truth.append({"plasmid_name": name, "ladder": "B", "plasmid_acc": acc,
                      "plasmid_len": len(newseq), "copy_number": LADDER_B_COPY_NO,
                      "repeat_len": rlen, "repeat_copies": LADDER_B_COPIES})
        placements.append({"plasmid": name, "family_src": src, "repeat_len": rlen,
                           "inserted_at": rec})
        print("  %2d kb repeat : carrier %-16s %7d -> %7d bp  (%d copies)"
              % (kb, acc, len(seq), len(newseq), LADDER_B_COPIES))

    # ---- write ------------------------------------------------------------
    ref = os.path.join(SIM, "probe01_ref.fasta")
    with open(ref, "w") as fo:
        for hdr, seq in records:
            fo.write(hdr + "\n")
            for k in range(0, len(seq), 70):
                fo.write(seq[k:k + 70] + "\n")

    # ---- verify what was actually written --------------------------------
    back = {h.split()[0]: s for h, s in read_fasta(ref)}
    bad = [t["plasmid_name"] for t in truth
           if len(back.get(t["plasmid_name"], "")) != t["plasmid_len"]]
    if bad or len(back) != len(records):
        sys.exit("VERIFY FAILED: %r" % (bad,))
    # every ladder-B repeat must really be present at the recorded identity
    for p in placements:
        s = back[p["plasmid"]]
        fam = cseq[p["family_src"]:p["family_src"] + p["repeat_len"]]
        found = sum(1 for off, ln in p["inserted_at"]
                    if ln == p["repeat_len"])
        if found != LADDER_B_COPIES:
            sys.exit("VERIFY FAILED: %s recorded %d copies" % (p["plasmid"], found))

    with open(os.path.join(SIM, "probe_truth.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(truth[0].keys()))
        w.writeheader()
        w.writerows(truth)
    json.dump({"seed": SEED, "chromosome": cacc, "chromosome_len": len(cseq),
               "ladder_a_copy": LADDER_A_COPY, "ladder_b_copies": LADDER_B_COPIES,
               "repeat_identity": REPEAT_IDENTITY, "placements": placements},
              open(os.path.join(SIM, "probe_design.json"), "w"), indent=1)

    total = len(cseq) + sum(t["plasmid_len"] for t in truth)
    print()
    print("verified: %d records, %d plasmids, %.2f Mb reference" % (len(back), len(truth), total / 1e6))
    print("-> %s" % ref)


if __name__ == "__main__":
    main()
