"""
Phase 0 / step 2 -- build the simulation arm's reference material.

COHORT v3. v2 (50 isolates, 200 assemblies, 11.6 h of Flye, zero failures)
passed 2 of its 4 gates and failed the other two. Both failures are fixed here,
and one of them is fixed by DELETING a mechanism rather than tuning it:

  (i)  CELL SPREAD 5.5 pp against a 10 pp floor. v2's 'small' band (4-12 kb x
       copy 8-30) is three regimes glued together, and only the middle one
       flips on prep. work/probe3 crossed 27 plasmids over 7.6-11.9 kb against
       copy {8,15,25} across all four gate cells and measured where the flips
       actually are. 'small' is now 8,400-10,900 (26.3 pp internal spread) and
       is 57 % of interventional slots instead of 44 %. Predicted gate spread
       15.0 pp -- stated here BEFORE the cohort is built.

  (ii) FRAGMENTED 19 against a floor of 20, with 0 of 148 induced. Three probes
       then showed the mechanism is not inducible in this simulator at all:
       repeat length 2-25 kb, 2-6 copies in plasmid and 2-4 in chromosome,
       identity 0.990 AND 1.000, r10_hq AND r9_raw, 30x AND 15x -> 0/164. Flye
       assembled every carrier complete and circular in all 16 cases. The
       'frag' role is WITHDRAWN and the floor with it. 'absorb' is untouched
       because it works (36/40 = 90 %).

KNOWN LIMITATION, RECORDED NOT HIDDEN. v2 held 45 real edges across 200 graphs,
175 of them edgeless, 91 % of nodes isolated -- and the three probes above,
built specifically to create branch points, produced one edge between them.
Badread samples uniformly from a perfect circular template and Flye traverses
that circle correctly whatever is inside it. This arm does not generate
assembly-graph topology, so it CANNOT test the relational-topology hypothesis
(B3/B4); that needs real assembly graphs. What it does support is pipeline
validation, the prep/depth contrast on absence, and 'absorbed'.

--------------------------------------------------------------------------
v2 was itself a rebuild of v1 (12 isolates), which exposed three problems:

  (a) SCALE. Leave-isolate-out CV and the paired bootstrap both resample
      ISOLATES, so v1's 168 label units were really n=12. B1-B0+ came out
      [+0.002,+0.136] in one arm and [-0.000,+0.142] in the other -- the same
      effect landing either side of the threshold on an arm that should not
      matter. N_ISOLATES is now 50 -> 200 assemblies.

  (b) STRATA. v1's 'small' band was 1-10 kb, which mixed two different failure
      modes. Sub-4 kb plasmids are shorter than nearly every read, so Flye reads
      the wrap-around as a repeat and drops them in ALL FOUR CELLS regardless of
      prep or depth -- those units are constants, not observations, and they
      dominated the positive class. They are now a separate 'tiny' band, tagged
      in ground_truth.csv, excluded from the interventional analysis and
      reported apart. 'small' is recentred on 4-12 kb, where dropout is
      probabilistic and where the ligation depletion actually operates.

  (c) REPEAT STRUCTURE. v1 realised recovered 119 / absent 32 / fragmented 8 /
      absorbed 1: the four-class mechanistic decomposition did not materialise.
      Fragmentation is repeat-mediated and absorption needs homology to a longer
      contig, but v1 drew plasmids independently from a pool so they shared no
      sequence with each other or with the chromosome -- nothing drove either
      mechanism. ~40% of isolates now carry an injected IS-like repeat family
      (see 'Repeat structure' below). This is the one axis where CONSTRUCTED
      isolates give leverage a real cohort cannot: the mechanisms can be induced
      on demand, and every placement is recorded in designs.json.

Three products:
  A) chromosome backbones   : the chromosome record of each ecoli_10 genome.
     NOTE: 6 of the 10 NCBI genomes are single-contig (chromosome only). That
     does NOT make them unusable. We are CONSTRUCTING isolates, so the genome's
     own plasmid complement is discarded either way -- we only ever take the
     chromosome and supply plasmids ourselves from the filtered pool. All 10
     are therefore usable backbones (each is reused across 5 isolates, with its
     own independent repeat injection, so no two isolates share a sequence).
  B) plasmid pool           : size-stratified draw from the filtered manifest,
     restricted to Enterobacterales hosts so the combination is plausible on an
     E. coli chromosome.
  C) per-isolate references : chromosome + k plasmids, with Badread
     'depth=' / 'circular=' headers, plus a ground-truth table.

Repeat structure
----------------
A repeat isolate gets one IS-like FAMILY: a 1-2 kb segment lifted from its own
chromosome (so base composition is genuinely genomic, and the family has a real
native locus). Copies are re-inserted at ~99% identity -- substitutions only, so
every copy is the same length.

v3 HAS ONE ROLE. The 'frag' role existed in v2 and is withdrawn; see (ii) above
for the 0/164 evidence. Its code path is gone from plan_repeats' caller, not
merely unused, so it cannot be re-enabled by accident.

  role 'absorb'  IS-mediated cointegrate: a ~99%-identity copy of the whole
                 plasmid is embedded in the chromosome, flanked by two copies of
                 the IS family, and the free plasmid is dropped to copy number 1.
                 Flye collapses the two near-identical copies, so >=95% of the
                 plasmid lands inside a contig >=3x its length -> 'absorbed'.

  DEVIATION FROM SPEC, STATED NOT BURIED: the brief asked for 1-2 kb shared
  segments only. Dispersed 1-2 kb elements cannot produce 'absorbed', because
  that rule needs >=95% of the PLASMID inside one contig, and only the shared
  bases align -- a 2 kb element in a 5 kb plasmid is 40% coverage, not 95%.
  Reaching 95% with 1-2 kb elements would require a plasmid that is a tandem
  array of ~20 copies, which is not IS-like biology. The cointegrate is the
  mechanism that actually produces absorption in nature and it is IS-mediated,
  so it is implemented as the 'absorb' role and flagged per placement
  (kind='cointegrate') in designs.json.

Ground truth is known by construction: we know exactly which plasmids went in,
at what copy number, with which repeat placements. That is what makes this arm
able to validate the labeller.

Reproducibility: everything derives from SEED. The pool draw uses one global
RNG; each isolate then uses its own RNG seeded SEED + ISOLATE_SEED_STRIDE + i,
so the isolate stream is independent of how many isolates precede it. Per-isolate
seeds, the composition, the repeat plan and every placement offset are written to
designs.json; ground_truth.csv carries the per-plasmid view of the same facts.
"""
import os, csv, random, collections, json

RAW       = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\Dataset\sequences.fasta"
MANIFEST  = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\work\refs\plasmid_manifest.csv"
ECOLI_DIR = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\Dataset\ecoli_10\ncbi_dataset\data"
REFS      = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\work\refs"
SIM       = r"D:\Plasmid-GNN\DNA_Sequencing_Technology\work\sim"

COHORT_VERSION      = 3
SEED                = 20260811
ISOLATE_SEED_STRIDE = 100_000      # isolate i uses SEED + STRIDE + i
POOL_PER_STRATUM    = 200
N_ISOLATES          = 50           # 200 assemblies; see (a) below

# Hosts plausible as donors for plasmids carried on an E. coli chromosome.
ENTERO = {"Escherichia", "Klebsiella", "Salmonella", "Enterobacter",
          "Citrobacter", "Serratia", "Shigella"}

# ---------------------------------------------------------------- strata
# COHORT v3 BOUNDS ARE MEASURED, NOT CHOSEN.
#
# v2 set 'small' to 4-12 kb x copy 8-30 and scored a 5.5 pp cell spread against
# a 10 pp floor, because that box is three regimes glued together: units that
# fail in BOTH preps, units that pass in BOTH, and the band that actually flips.
# Only a unit that flips contributes to the spread.
#
# probe3 (work/probe3) crossed 27 plasmids over 7.6-11.9 kb against copy
# {8, 15, 25} and assembled all four gate cells. Failure rate per cell, and the
# resulting internal spread, by candidate band:
#
#     band            n     L30   R30   L15   R15    spread
#     1,000- 8,400     5   80.0  60.0  60.0  60.0    20.0 pp   fails in both preps
#     7,500-12,500    27   33.3  18.5  33.3  14.8    18.5 pp   the whole probe
#     8,400-10,900    19   26.3  10.5  31.6   5.3    26.3 pp   <- chosen
#     9,000-10,500    13   15.4   0.0  23.1   0.0    23.1 pp   n too small
#
# 8,400-10,900 is chosen over the narrower windows because its 26.3 pp rests on
# 19 plasmids rather than a handful. The report script's own suggestion
# (8,828-9,251, "100 % flip rate") is min/max over the 3 plasmids that flipped
# at both depths and is overfitted -- it is not used.
#
# 'tiny' now runs to 8,400, so it is no longer literally tiny: it is THE
# DETERMINISTIC BAND, everything that fails in both preps at this read length.
# The key name is kept because ground_truth.csv, 04b and 09_report.py all index
# on it and a rename would ripple through four files for no analytic gain.
#
# These bounds are valid for FRAGLEN 15000,13000 ONLY. Change the fragment
# length and every one of them moves -- re-run probe_v3c.* before trusting them.
STRATA = ("tiny", "small", "medium", "large")
STRATUM_BOUNDS = {
    "tiny":   (  1_000,   8_400),   # deterministic band -- excluded, reported apart
    "small":  (  8_400,  10_900),   # THE prep window: 26.3 pp internal spread
    "medium": ( 10_900, 100_000),
    "large":  (100_000, 500_000),
}
ANALYSIS_SET = {
    "tiny":   "deterministic_tiny",   # excluded from interventional analysis
    "small":  "interventional",
    "medium": "interventional",
    "large":  "interventional",
}

# Copy number by stratum: short plasmids are high-copy, large ones low-copy.
# 'small' is capped at 25 rather than 30 because probe3 crossed copy {8, 15, 25}
# and flips occurred at ALL THREE levels (1 / 1 / 1 of the both-depth flips), so
# copy number is not the discriminator inside the band -- but 30 was never
# tested and the range is kept inside what was measured.
COPY_RANGE = {"tiny": (10, 40), "small": (8, 25), "medium": (2, 6), "large": (1, 2)}

# ---------------------------------------------------------------- repeats
REPEAT_ISOLATE_FRACTION = 0.40     # ~40% of isolates carry an IS family
IS_LEN_RANGE            = (1_000, 2_000)
IS_IDENTITY             = 0.99     # substitutions only -> copies keep their length
FRAG_CHROM_COPIES       = (2, 4)
FRAG_PRIMARY_COPIES     = (2, 3)   # copies inside the main carrier plasmid
FRAG_SECONDARY_COPIES   = (1, 2)   # copies inside a second plasmid (plasmid<->plasmid)
ABSORB_CHROM_COPIES     = (2, 3)
COINTEGRATE_COPY_NUMBER = 1        # free plasmid drops to 1x when cointegrated

# Composition sweep: 25 patterns, cycled twice over 50 isolates, so plasmid
# count and stratum mix move independently of each other and of the chemistry
# arm.
#
# v3 REBALANCE. Only units inside the prep window contribute to the cell-spread
# gate; medium and large fail ~6 % in every cell and contribute ~0. v2 drew
# small at 32 of 73 interventional slots = 44 %, so its 18.5 pp of internal
# spread arrived at the gate as 8.1 pp. 'small' is now the majority stratum.
#
# Per cycle: tiny 6 / small 45 / medium 23 / large 11 = 85 plasmids, so 170 over
# the cohort and 158 interventional. Small is 45 / 79 = 57 % of interventional
# slots, and 0.57 x 26.3 pp = 15.0 pp expected at the gate against a 10 pp floor.
#
# 'tiny' is cut from 14 to 6 per cycle: it is excluded from the interventional
# analysis by construction, so slots spent there buy nothing but a bigger
# reported-apart table. Enough is kept to keep reporting it honest.
COMPOSITIONS = [
    ("small", "small", "medium"),
    ("small", "small", "medium", "large"),
    ("tiny", "small", "small", "medium"),
    ("small", "small", "large"),
    ("small", "medium"),
    ("small", "small", "medium", "medium"),
    ("tiny", "small", "small", "large"),
    ("small", "small", "medium"),
    ("small", "medium", "large"),
    ("small", "small", "medium", "large"),
    ("tiny", "small", "small", "medium"),
    ("small", "small", "medium"),
    ("small", "small", "large"),
    ("small", "medium", "medium"),
    ("tiny", "small", "small", "medium"),
    ("small", "small", "medium", "large"),
    ("small", "small", "medium"),
    ("small", "medium"),
    ("tiny", "small", "small", "large"),
    ("small", "small", "medium", "medium"),
    ("small", "small", "medium"),
    ("small", "medium", "large"),
    ("tiny", "small", "small", "medium"),
    ("small", "small", "large"),
    ("small", "small", "medium", "large"),
]

BASES = b"ACGT"


def stratum(L):
    for s in STRATA:
        lo, hi = STRATUM_BOUNDS[s]
        if lo <= L < hi:
            return s
    return "large" if L >= STRATUM_BOUNDS["large"][0] else "tiny"


def read_fasta(path):
    hdr, chunks = None, []
    with open(path, "rb") as fh:
        for line in fh:
            if line[:1] == b">":
                if hdr is not None:
                    yield hdr, b"".join(chunks)
                hdr = line[1:].decode("utf-8", "replace").rstrip()
                chunks = []
            else:
                chunks.append(line.strip())
    if hdr is not None:
        yield hdr, b"".join(chunks)


def wrap(seq, width=70):
    return b"\n".join(seq[i:i + width] for i in range(0, len(seq), width))


# ---------------------------------------------------------------- repeat helpers
def mutate(seq, rng, identity):
    """Substitution-only copy at approximately `identity`.

    Substitutions only (no indels) so every copy is exactly len(seq) -- that
    keeps the offset bookkeeping below exact, and it is what a recently
    transposed IS copy looks like anyway.
    """
    b = bytearray(seq)
    n_sub = int(round(len(b) * (1.0 - identity)))
    if n_sub and len(b):
        for pos in rng.sample(range(len(b)), min(n_sub, len(b))):
            cur = b[pos]
            alt = [c for c in BASES if c != cur]
            b[pos] = rng.choice(alt)
    realised = sum(1 for i in range(len(b)) if b[i] == seq[i]) / len(b) if len(b) else 1.0
    return bytes(b), realised


def insert_copies(seq, items, replicon, family_id, out):
    """Insert every item into `seq` in ONE pass and record exact final offsets.

    items = [(position, payload_bytes, identity, kind, extra_dict)], positions
    in ORIGINAL coordinates and distinct.

    Inserted descending so an insertion never shifts a position not yet used;
    the FINAL offset of each copy is then its original position plus the total
    length inserted strictly before it. This has to be one pass per replicon: a
    second round of insertions would silently invalidate the offsets recorded by
    the first, which is exactly the bug this signature exists to prevent. Both
    offsets are recorded so designs.json is auditable against the written FASTA.
    """
    for pos, payload, _i, _k, _x in sorted(items, key=lambda t: t[0], reverse=True):
        seq = seq[:pos] + payload + seq[pos:]
    for ci, (pos, payload, ident, kind, extra) in enumerate(
            sorted(items, key=lambda t: t[0])):
        shift = sum(len(p) for q, p, _i, _k, _x in items if q < pos)
        rec = {"replicon": replicon, "family": family_id, "kind": kind,
               "copy_index": ci, "length": len(payload),
               "identity": round(ident, 5),
               "offset_original": pos, "offset_final": pos + shift}
        rec.update(extra)
        out.append(rec)
    return seq


def pick_positions(rng, seq_len, n, margin=200):
    """n distinct insertion points, kept away from the sequence ends."""
    lo, hi = margin, max(margin + 1, seq_len - margin)
    if hi - lo <= n:
        return sorted(rng.randrange(0, max(1, seq_len)) for _ in range(n))
    return sorted(rng.sample(range(lo, hi), n))


# ---------------------------------------------------------------- products A/B
def build_chromosomes():
    """Largest record of each genome = chromosome."""
    os.makedirs(REFS, exist_ok=True)
    out = {}
    for acc in sorted(os.listdir(ECOLI_DIR)):
        d = os.path.join(ECOLI_DIR, acc)
        if not os.path.isdir(d):
            continue
        fna = [f for f in os.listdir(d) if f.endswith(".fna")]
        if not fna:
            continue
        recs = list(read_fasta(os.path.join(d, fna[0])))
        hdr, seq = max(recs, key=lambda r: len(r[1]))
        if len(seq) < 500_000:
            print("  SKIP %s: largest record %d bp < 500 kb" % (acc, len(seq)))
            continue
        out[acc] = (hdr, seq)
        print("  %s  chromosome %9d bp  (%d records in file, %d discarded)"
              % (acc, len(seq), len(recs), len(recs) - 1))
    return out


def select_pool(rng):
    rows = [r for r in csv.DictReader(open(MANIFEST, encoding="utf-8"))
            if r["verdict"] == "keep" and r["genus"] in ENTERO]
    by = collections.defaultdict(list)
    for r in rows:
        by[stratum(int(r["length"]))].append(r)
    print("  candidates (Enterobacterales, post-filter):")
    for s in STRATA:
        lo, hi = STRATUM_BOUNDS[s]
        print("     %-7s %-18s %d" % (s, "[%d, %d)" % (lo, hi), len(by[s])))
    chosen = []
    for s in STRATA:
        pool = by[s]
        if len(pool) < POOL_PER_STRATUM:
            print("     !! only %d candidates in '%s', wanted %d"
                  % (len(pool), s, POOL_PER_STRATUM))
        rng.shuffle(pool)
        chosen.extend(pool[:POOL_PER_STRATUM])
    return {r["accession"]: r for r in chosen}


def extract_pool_seqs(wanted):
    """One streaming pass over the 6.9 GB raw file for the chosen accessions."""
    got = {}
    for hdr, seq in read_fasta(RAW):
        acc = hdr.split()[0]
        if acc in wanted and acc not in got:
            got[acc] = seq
            if len(got) == len(wanted):
                break
    return got


# ---------------------------------------------------------------- product C
def plan_repeats(role, rng, comp):
    """Decide this isolate's repeat plan given its assigned role.

    The plan is fixed BEFORE plasmids are drawn so the draw can respect each
    slot's insertion budget and keep the plasmid inside its stratum.
    """
    fam_len = rng.randint(*IS_LEN_RANGE)
    plan = {"role": role, "family_len": fam_len, "identity_target": IS_IDENTITY,
            "slot_copies": {}, "cointegrate_slot": None}

    if role == "frag":
        plan["chrom_copies"] = rng.randint(*FRAG_CHROM_COPIES)
        # carrier = a mid-size slot if one exists, else the largest available
        prefer = [i for i, s in enumerate(comp) if s in ("small", "medium")]
        if not prefer:
            prefer = list(range(len(comp)))
        primary = rng.choice(prefer)
        plan["slot_copies"][primary] = rng.randint(*FRAG_PRIMARY_COPIES)
        others = [i for i in range(len(comp)) if i != primary]
        if others:
            secondary = rng.choice(others)
            plan["slot_copies"][secondary] = rng.randint(*FRAG_SECONDARY_COPIES)
    else:
        plan["chrom_copies"] = rng.randint(*ABSORB_CHROM_COPIES)
        # v3: cointegrate a MEDIUM slot by preference, not the shortest.
        #
        # A cointegrated plasmid is `absorbed` in all four cells -- it is a
        # constant, and constants contribute nothing to the cell-spread gate.
        # v2 cointegrated the shortest slot, which was always 'small', so 40 of
        # the prep-sensitive units were spent on a class that cannot flip.
        # Taking them from 'medium' instead leaves the prep window intact.
        #
        # Medium satisfies the absorbed rule just as well: the test is >=95 % of
        # the plasmid inside one contig that is >=3x its length, and a 10.9-100 kb
        # plasmid embedded in a ~4.8 Mb chromosome clears 3x with room to spare.
        order = {"medium": 0, "large": 1, "small": 2, "tiny": 3}
        cand = sorted(range(len(comp)), key=lambda i: (order[comp[i]], i))
        plan["cointegrate_slot"] = cand[0]
    return plan


def draw_plasmid(avail, seqs, st, budget):
    """Pop an accession from stratum `st` that stays in band after `budget`
    inserted bases. Scans from the tail so budget=0 reproduces a plain pop()."""
    lo, hi = STRATUM_BOUNDS[st]
    lst = avail[st]
    for k in range(len(lst) - 1, -1, -1):
        if lo <= len(seqs[lst[k]]) + budget < hi:
            return lst.pop(k), True
    return (lst.pop() if lst else None), False


def build_isolate(i, iso, arm, ck, chdr, cseq, comp, plan, avail, seqs, wanted):
    """Materialise one isolate. Returns (records, truth_rows, design)."""
    rng = random.Random(SEED + ISOLATE_SEED_STRIDE + i)
    placements = []
    family_id = family_seq = None
    fam_src = None

    if plan is not None:
        fam_len = plan["family_len"]
        fam_src = rng.randrange(0, max(1, len(cseq) - fam_len))
        family_seq = cseq[fam_src:fam_src + fam_len]
        family_id = "IS_%s" % iso

    # --- draw plasmids, respecting each slot's insertion budget
    picks, out_of_band = [], 0
    for slot, st in enumerate(comp):
        budget = 0
        if plan is not None:
            budget = plan["slot_copies"].get(slot, 0) * plan["family_len"]
        acc, in_band = draw_plasmid(avail, seqs, st, budget)
        if acc is None:
            continue
        if not in_band:
            out_of_band += 1
        picks.append((slot, acc))

    # --- inject repeat copies into plasmids
    plasmid_seqs, repeat_meta = {}, {}
    for slot, acc in picks:
        s = seqs[acc]
        n_copies = plan["slot_copies"].get(slot, 0) if plan is not None else 0
        if n_copies:
            pos = pick_positions(rng, len(s), n_copies)
            items = []
            for p in pos:
                cp, idt = mutate(family_seq, rng, IS_IDENTITY)
                items.append((p, cp, idt, "is_copy", {}))
            s = insert_copies(s, items, "plasmid_slot_%d" % slot,
                              family_id, placements)
        plasmid_seqs[slot] = s
        repeat_meta[slot] = {"copies": n_copies,
                             "bases": n_copies * (plan["family_len"] if plan else 0)}

    # --- inject into the chromosome (ONE pass: see insert_copies docstring)
    cointegrated_slot = None
    if plan is not None:
        cs = plan["cointegrate_slot"]
        if cs is not None and cs in plasmid_seqs:
            cointegrated_slot = cs
        n_pos = plan["chrom_copies"] + (1 if cointegrated_slot is not None else 0)
        pos = pick_positions(rng, len(cseq), n_pos)
        items = []
        for p in pos[:plan["chrom_copies"]]:
            cp, idt = mutate(family_seq, rng, IS_IDENTITY)
            items.append((p, cp, idt, "is_copy", {}))
        if cointegrated_slot is not None:
            # IS-mediated cointegrate: [IS][~99% copy of the plasmid][IS]
            left, _ = mutate(family_seq, rng, IS_IDENTITY)
            right, _ = mutate(family_seq, rng, IS_IDENTITY)
            body, body_idt = mutate(plasmid_seqs[cointegrated_slot], rng, IS_IDENTITY)
            items.append((pos[-1], left + body + right, body_idt, "cointegrate",
                          {"cointegrated_slot": cointegrated_slot,
                           "plasmid_body_length": len(body),
                           "flank_length": len(left),
                           "identity_is": "identity refers to the plasmid body"}))
        cseq = insert_copies(cseq, items, "chromosome", family_id, placements)

    # --- write the reference FASTA
    path = os.path.join(SIM, iso + "_ref.fasta")
    truth_rows = []
    with open(path, "wb") as fo:
        fo.write(b">chromosome depth=1.0 circular=true\n")
        fo.write(wrap(cseq) + b"\n")
        for j, (slot, acc) in enumerate(picks, 1):
            s = plasmid_seqs[slot]
            st = stratum(len(s))
            cointegrated = cointegrated_slot == slot
            if cointegrated:
                cn = COINTEGRATE_COPY_NUMBER
            else:
                lo, hi = COPY_RANGE[st]
                cn = rng.randint(lo, hi)
            name = "plasmid_%d" % j
            fo.write((">%s depth=%d circular=true\n" % (name, cn)).encode())
            fo.write(wrap(s) + b"\n")
            rb = repeat_meta[slot]["bases"]
            truth_rows.append({
                "isolate": iso, "arm": arm, "chromosome_acc": ck,
                "chromosome_len": len(cseq), "plasmid_name": name,
                "plasmid_acc": acc, "plasmid_len": len(s),
                "plasmid_len_original": len(seqs[acc]),
                "stratum": st, "stratum_designed": comp[slot],
                "analysis_set": ANALYSIS_SET[st], "copy_number": cn,
                "organism": wanted[acc]["organism"], "gc": wanted[acc]["gc"],
                "repeat_family": family_id or "",
                "repeat_role": (plan["role"] if plan else ""),
                "repeat_copies_in_plasmid": repeat_meta[slot]["copies"],
                "repeat_bases_in_plasmid": rb,
                "repeat_frac": round(rb / len(s), 5) if len(s) else 0.0,
                "cointegrated": "true" if cointegrated else "false",
                "isolate_seed": SEED + ISOLATE_SEED_STRIDE + i,
            })

    design = {
        "isolate": iso, "arm": arm, "chromosome": ck,
        "chromosome_len": len(cseq), "n_plasmids": len(picks),
        "composition": list(comp), "ref": path,
        "isolate_seed": SEED + ISOLATE_SEED_STRIDE + i,
        "plasmid_accessions": [a for _, a in picks],
        "slots_out_of_band": out_of_band,
        "repeat_structure": None if plan is None else {
            "family_id": family_id, "family_length": plan["family_len"],
            "family_source_offset_in_backbone": fam_src,
            "identity_target": IS_IDENTITY, "role": plan["role"],
            "chromosome_copies": plan["chrom_copies"],
            "slot_copies": {str(k): v for k, v in plan["slot_copies"].items()},
            "cointegrate_slot": cointegrated_slot,
            "placements": placements,
        },
    }
    return truth_rows, design


def check_stale_cohort():
    """Refuse to overwrite the references while a previous cohort's assemblies
    are still on disk.

    This is the one way the rebuild fails silently instead of loudly. Nothing
    downstream is keyed to a cohort version: 04_label.py walks every directory
    in asm/ and labels it against sim/<isolate>_ref.fasta, and 05_features.py
    walks the same list. After this script runs, sim01_ref.fasta is a DIFFERENT
    ISOLATE with different plasmids -- so a surviving v1 asm/sim01_lig_d15 would
    be labelled against plasmids that were never in it, come back all-absent,
    and flow into node_features.csv and the ladder as if it were real. The
    accounting check in 04b would flag it, but only after everything downstream
    had already consumed it.

    Set REBUILD_FORCE=1 to proceed anyway (only sensible if asm/ and graphs/ are
    already empty and this is a re-run of this script alone).
    """
    root = os.path.dirname(SIM)
    stale = []
    for sub in ("asm", "graphs"):
        d = os.path.join(root, sub)
        if os.path.isdir(d):
            n = len([x for x in os.listdir(d)
                     if os.path.isdir(os.path.join(d, x))])
            if n:
                stale.append((sub, n))
    old_refs = []
    if os.path.isdir(SIM):
        for f in sorted(os.listdir(SIM)):
            if f.endswith("_ref.fasta"):
                old_refs.append(f)
    if not stale and not old_refs:
        return
    if os.environ.get("REBUILD_FORCE") == "1":
        print("  REBUILD_FORCE=1 -- proceeding over %d stale ref(s) and %s"
              % (len(old_refs), stale or "no assemblies"))
        return
    print("=" * 66)
    print("REFUSING TO REBUILD -- a previous cohort is still on disk")
    print("=" * 66)
    for sub, n in stale:
        print("  %-8s %d directories" % (sub + "/", n))
    print("  sim/     %d existing *_ref.fasta" % len(old_refs))
    print()
    print("  Every isolate's plasmid complement changes in this rebuild, so a")
    print("  surviving assembly would be labelled against plasmids it never")
    print("  contained. Clear the old cohort first:")
    print()
    print("     rm -rf %s %s" % (os.path.join(root, "asm"),
                                 os.path.join(root, "graphs")))
    print("     rm -f  %s" % os.path.join(SIM, "sim*_ref.fasta"))
    print("     rm -f  %s" % os.path.join(root, "results", "node_features.csv"))
    print("     rm -f  %s" % os.path.join(root, "results", "graph_diagnostics.json"))
    print("     rm -f  %s" % os.path.join(root, "results", "pyg_dataset.pt"))
    print("     rm -f  %s" % os.path.join(root, "results", "label_validation.csv"))
    print("     bash scripts/reset_scratch.sh     # drop any half-built read sets")
    print()
    print("  Then rerun this script. (REBUILD_FORCE=1 overrides this check.)")
    raise SystemExit(2)


def main():
    check_stale_cohort()
    rng = random.Random(SEED)
    os.makedirs(SIM, exist_ok=True)

    print("=" * 66)
    print("COHORT v%d -- %d isolates x {lig,rap} x {30x,15x} = %d assemblies"
          % (COHORT_VERSION, N_ISOLATES, N_ISOLATES * 4))
    print("=" * 66)

    print()
    print("A) chromosome backbones")
    print("-" * 66)
    chroms = build_chromosomes()
    print("  -> %d usable backbones (each reused across %d isolates, with"
          % (len(chroms), N_ISOLATES // max(1, len(chroms))))
    print("     independent repeat injection, so no two references are identical)")

    print()
    print("B) plasmid pool")
    print("-" * 66)
    wanted = select_pool(rng)
    print("  selected %d accessions; extracting sequences..." % len(wanted))
    seqs = extract_pool_seqs(wanted)
    print("  extracted %d / %d" % (len(seqs), len(wanted)))

    pool_fa = os.path.join(REFS, "plasmid_pool.fasta")
    with open(pool_fa, "wb") as fo:
        for acc, s in seqs.items():
            m = wanted[acc]
            fo.write((">%s %s len=%d\n" % (acc, m["organism"], len(s))).encode())
            fo.write(wrap(s) + b"\n")
    print("  -> %s" % pool_fa)

    print()
    print("C) per-isolate references")
    print("-" * 66)
    avail = collections.defaultdict(list)
    for acc, s in seqs.items():
        avail[stratum(len(s))].append(acc)
    for s in avail:
        rng.shuffle(avail[s])

    # required draw per stratum, so a shortfall is caught before any file is written
    need = collections.Counter()
    for i in range(N_ISOLATES):
        for st in COMPOSITIONS[i % len(COMPOSITIONS)]:
            need[st] += 1
    short = [s for s in STRATA if len(avail[s]) < need[s]]
    print("  draw required / available:")
    for s in STRATA:
        print("     %-7s %3d / %3d%s"
              % (s, need[s], len(avail[s]), "   SHORTFALL" if s in short else ""))
    if short:
        raise SystemExit("  !! pool too small for strata %s -- raise "
                         "POOL_PER_STRATUM and rerun" % short)

    # Repeat-carrying isolates, balanced across chemistry arms so the repeat
    # axis is not confounded with chemistry.
    #
    # COHORT v3: EVERY repeat isolate now takes the 'absorb' role. The 'frag'
    # role is withdrawn -- it is not a tuning change, the mechanism does not
    # exist in this simulator. Measured:
    #
    #   cohort v2      IS carrier -> fragmented   0 / 148
    #   probe round 1  repeats of 3, 8, 15, 25 kb at 99 % identity, 2 copies,
    #                  carriers 94-200 kb                          0 / 4
    #   probe round 2  repeats of 2-25 kb at 100 % identity (verified verbatim
    #                  in the written FASTA), 2-6 copies in the plasmid and
    #                  2-4 in the chromosome, under r10_hq @ 30x AND
    #                  r9_raw @ 15x                                0 / 12
    #
    # Flye assembled every carrier complete and circular in all 16 cases, and
    # emitted the repeats as separate high-coverage contigs beside them -- it
    # detects them and resolves them anyway. Badread samples uniformly from a
    # perfect circular template; real fragmentation comes from coverage
    # dropouts, chimeras and strain heterogeneity it does not reproduce.
    # Depth is not the lever either: in v2, medium/large units below 1.5x
    # relative depth came out 56 recovered vs 6 fragmented.
    #
    # 'absorb' is kept unchanged because it WORKS: 36/40 = 90 % conversion.
    # Absorption is a containment relationship and needs no read to span
    # anything, which is exactly why it survives where 'frag' does not.
    n_rep = int(round(REPEAT_ISOLATE_FRACTION * N_ISOLATES))
    even = [i for i in range(N_ISOLATES) if i % 2 == 0]     # r10_hq
    odd = [i for i in range(N_ISOLATES) if i % 2 == 1]      # r9_raw
    rng.shuffle(even); rng.shuffle(odd)
    role_of = {}
    for group, k in ((even, n_rep // 2), (odd, n_rep - n_rep // 2)):
        for ix in sorted(group[:k]):
            role_of[ix] = "absorb"
    rep_idx = sorted(role_of)

    chrom_keys = sorted(chroms)
    truth_rows, designs = [], []

    for i in range(N_ISOLATES):
        iso = "sim%02d" % (i + 1)
        ck = chrom_keys[i % len(chrom_keys)]
        chdr, cseq = chroms[ck]
        comp = COMPOSITIONS[i % len(COMPOSITIONS)]
        arm = "r10_hq" if i % 2 == 0 else "r9_raw"

        plan = None
        if i in role_of:
            # planning uses its own derived stream so it cannot consume from the
            # stream build_isolate() replays for the sequence construction
            plan = plan_repeats(role_of[i],
                                random.Random(SEED + 2 * ISOLATE_SEED_STRIDE + i),
                                comp)

        rows, design = build_isolate(i, iso, arm, ck, chdr, cseq, comp, plan,
                                     avail, seqs, wanted)
        truth_rows.extend(rows)
        designs.append(design)
        rs = design["repeat_structure"]
        print("  %s arm=%-7s chrom=%s %d plasmids (%s)%s"
              % (iso, arm, ck, design["n_plasmids"], ",".join(comp),
                 "  repeats=%s/%dbp" % (rs["role"], rs["family_length"]) if rs else ""))

    # ------------------------------------------------------------ manifests
    tp = os.path.join(SIM, "ground_truth.csv")
    with open(tp, "w", newline="", encoding="utf-8") as fo:
        w = csv.DictWriter(fo, fieldnames=list(truth_rows[0].keys()))
        w.writeheader()
        w.writerows(truth_rows)

    provenance = {
        "cohort_version": COHORT_VERSION,
        "seed": SEED,
        "isolate_seed_stride": ISOLATE_SEED_STRIDE,
        "n_isolates": N_ISOLATES,
        "pool_per_stratum": POOL_PER_STRATUM,
        "strata": {s: list(STRATUM_BOUNDS[s]) for s in STRATA},
        "analysis_set": ANALYSIS_SET,
        "copy_range": {s: list(COPY_RANGE[s]) for s in STRATA},
        "repeat": {
            "isolate_fraction": REPEAT_ISOLATE_FRACTION,
            "n_repeat_isolates": len(rep_idx),
            "is_len_range": list(IS_LEN_RANGE),
            "is_identity": IS_IDENTITY,
            "cointegrate_copy_number": COINTEGRATE_COPY_NUMBER,
            "repeat_isolates": ["sim%02d" % (i + 1) for i in rep_idx],
        },
        "compositions": [list(c) for c in COMPOSITIONS],
        "inputs": {"raw_plasmids": RAW, "manifest": MANIFEST,
                   "chromosome_dir": ECOLI_DIR},
        "cells_per_isolate": ["lig_d30", "lig_d15", "rap_d30", "rap_d15"],
    }
    dp = os.path.join(SIM, "designs.json")
    json.dump({"provenance": provenance, "isolates": designs},
              open(dp, "w"), indent=2)

    # ------------------------------------------------------------ summary
    by_stratum = collections.Counter(r["stratum"] for r in truth_rows)
    by_set = collections.Counter(r["analysis_set"] for r in truth_rows)
    roles = collections.Counter(d["repeat_structure"]["role"] for d in designs
                                if d["repeat_structure"])
    n_coint = sum(1 for r in truth_rows if r["cointegrated"] == "true")
    n_carrier = sum(1 for r in truth_rows if int(r["repeat_copies_in_plasmid"]) > 0)
    oob = sum(d["slots_out_of_band"] for d in designs)

    print()
    print("=" * 66)
    print("COHORT SUMMARY")
    print("=" * 66)
    print("  isolates                 : %d  (%d assemblies at 4 cells each)"
          % (len(designs), 4 * len(designs)))
    print("  plasmid label units      : %d  (x4 cells = %d)"
          % (len(truth_rows), 4 * len(truth_rows)))
    print("  strata (realised)        :", {s: by_stratum[s] for s in STRATA})
    print("  analysis sets            :", dict(by_set))
    print("  interventional units     : %d  (x4 cells = %d)"
          % (by_set["interventional"], 4 * by_set["interventional"]))
    print("  arms                     :",
          dict(collections.Counter(d["arm"] for d in designs)))
    print("  repeat isolates          : %d  roles=%s" % (len(rep_idx), dict(roles)))
    print("  plasmids carrying IS     : %d" % n_carrier)
    print("  cointegrated plasmids    : %d   <- the 'absorbed' candidates" % n_coint)
    print("  slots drawn out of band  : %d" % oob)
    print()
    print("  EXPECTED-CLASS CHECK (floors are enforced later, in 09_report.py):")
    print("    absorbed floor is 20 and there are %d cointegrates -- one per"
          % n_coint)
    print("    absorb-role isolate, each seen in all 4 cells. v2 converted 90%.")
    print("    fragmented has NO floor from v3 on. The 'frag' role is withdrawn:")
    print("    0/164 induced across the v2 cohort and three probes. Whatever")
    print("    fragmented units appear are stochastic background (v2: ~3%) and")
    print("    carry no claim -- see 09_report.py CLASS_FLOORS.")
    print()
    print("  -> %s" % tp)
    print("  -> %s" % dp)


if __name__ == "__main__":
    main()
