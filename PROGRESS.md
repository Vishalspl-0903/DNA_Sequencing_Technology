# Plasmid Recovery GNN — Phase 0 progress and results

**Date:** 2026-08-12
**Scope of this document:** what has been built, what was verified, what the
numbers actually say, and what is still running or not yet started.

Status key: **[DONE]** verified · **[RUNNING]** in progress · **[PENDING]** not started
· **[SUPERSEDED]** correct when written, replaced by a later section

---

## 0. One-paragraph summary

The proposal's factual claims about the Johnson cohort all check out against the
deposited data. The dataset on disk, however, was label/metadata scaffolding
with **no assembly graphs at all** — the model's entire input modality was
missing. Priority therefore moved from acquiring more real data to standing up
the simulation arm end to end. A working toolchain, a filtered plasmid reference
set, a simulated cohort, a unit-tested labeller, GFA feature extraction and PyG
graph export are in place and validated.

**The first full cohort ran, and it failed its own adequacy checks.** The
pipeline was not at fault — zero Flye failures, the labeller unit-tested, every
designed plasmid labelled exactly once. The **cohort design** was: 12 isolates
is n=12 for every test that matters, the 1–10 kb size band mixed a
probabilistic failure mode with a deterministic one, and plasmids drawn
independently from a pool share no sequence, so the two repeat-mediated label
classes never materialised (fragmented 8, absorbed 1). Sec. 10 below specifies
the v2 rebuild that fixes all three, and Sec. 11 records four numbered gates
that must pass before the GNN is trained on anything.

---

## 1. Verification of the proposal against the deposited data **[DONE]**

### 1.1 Claims that check out

| Proposal claim | Data on disk | Verdict |
|---|---|---|
| 33 plasmids, 14 isolates, 1.9–194 kb | 33 / 14 / 1,919–194,062 bp | ✅ |
| Unicycler recovered 100 % | 99/99 rows `present` | ✅ |
| Flye recovers 88–91 % | flye-hq 88.9 %, flye-meta 90.9 %, flye-raw 90.9 % | ✅ |
| All 11 plasmids <10 kb missed by ≥1 long-read-only assembler | exactly 11 exist, all 11 missed | ✅ |
| Size bands <10 / 10–99 / ≥100 kb | 11 / 16 / 6 — maps exactly onto Johnson's small/medium/large | ✅ |
| Triplicate assemblies | 7 assemblers × 3 reps × 33 plasmids = 693 rows | ✅ |
| Pooled positive rate 9–12 % | flye-raw rep 1 = 3/33 = 9.1 % | ✅ |
| Wick = 7 isolates, deposited FASTA | `assemblies.tar.gz`, 7 references, 3–7 replicons each | ✅ |

The Wick references are well suited to the hypothesis: 6 of 7 carry a <10 kb
plasmid, one carries five.

### 1.2 Errors found

**a) The correctness gate cannot pass as written.** Under the pinned mode
(`--nano-raw` → Johnson's `flye-raw`), recovery is **30/33, 29/33, 31/33** across
replicates. A "±1 plasmid" tolerance around an unspecified target is *narrower
than the paper's own replicate spread*, so the gate passes or fails depending on
which replicate is silently chosen.

**b) Open item (iv) resolves negatively.** Genus distribution across the 92
assembled NEKSUS isolates:

| genus | n |
|---|---|
| *Escherichia* | 58 |
| *Klebsiella* | 29 |
| *Citrobacter* | 2 |
| *Enterobacter* | 2 |
| *Serratia* | 1 |

Only two genera clear 12 isolates, so a 12/12 genus-disjoint split can only ever
be E. coli vs Klebsiella. The remaining five isolates cannot form a side.

**c) Open item (i) resolves negatively.** Johnson's supplement contains **no
reference FASTAs** — two CSV tables only. Unicycler re-derivation is the
mainline, not a contingency, and it needs Johnson Illumina reads that the ~1.8 GB
transfer budget does not cover.

**d) NEKSUS is 96 deposited, 92 assembled.** The four without assemblies
(AF17, AF25, AF33, AHB9) include **three of the four** that the metadata flags
"likely contaminated with Escherichia coli" — an exclusion criterion the
proposal was inheriting silently.

**e) Table I mislabels the NEKSUS references.** What exists is
`contigs_summary_sup_cleaned.csv` (3,917 rows of per-contig metadata across six
assemblers), not FASTA. It also carries `circular` and `completeness` columns —
the same circularity leak Sec. IV excludes, in a more tempting form, since it is
the one file that invites direct feature extraction. **A parser-level guard now
exists** (`scripts/leak_guard.py`, Sec. 12.4) and it strips **five** columns from
that file, not two: `circular`, `completeness`, `no_circular_chromosomes`,
`no_circular_contigs` and `chromosome_circularity` are the same leak under
different names.

### 1.3 Problems inside the files that are present

- **`supplementary-data-1.csv` is partially corrupt**: 120 rows of which **60 are
  exact duplicates** (each Long Reads row triplicated with identical values — no
  usable replicate dimension). Covers 12 isolates / 30 plasmids;
  `2022NG-0076_1`, `2022QW-00133_2`, `2022QW-00133_3` are absent.
- **`sequences.fasta` is not PLSDB** — an unfiltered NCBI plasmid query: 72,556
  records, 7.3 Gbp, **1,581 records ≥500 kb**, largest **11.85 Mb**. At 6.92 GB it
  is ~23× the 300 MB budgeted for the whole simulation arm.
- **`ecoli_ids.txt`** holds 10 RefSeq accessions; `ecoli_10/` holds those genomes.

### 1.4 What was missing entirely

No reads (Wick's barcoded FASTQs, Johnson, NEKSUS), no Illumina, **no GFA files
at all**, no Johnson or NEKSUS reference FASTAs, no tooling, no code.

### 1.5 One encouraging pre-Phase-0 signal

From the NEKSUS contig summary alone, against `hybracter_hybrid` as reference
(353 plasmid-like contigs over 88 samples): `hybracter_long` is short by **39**
across 20 samples, `flye` short by **26** across 12 — and **all 26 flye deficits
fall in samples carrying a <10 kb plasmid, none outside**. A contig-count proxy,
not alignment-based labels, but it points the right way.

---

## 2. Toolchain **[DONE]**

Installed into WSL userspace via micromamba — **no sudo required** (the machine's
WSL user has no passwordless sudo, and the image lacks `bzip2`, so the
micromamba tarball is extracted with Python's `tarfile`).

| tool | version |
|---|---|
| Flye | 2.9.6-b1802 |
| minimap2 | 2.31-r1302 |
| seqkit | 2.13.0 |
| rasusa | 5.1.0 |
| Badread | 0.4.2 |
| PyTorch / PyG (Windows) | 2.7.1+cu118 / 2.8.0 |

**Corrections to the earlier version of this section**, both found while auditing
the working directory on 2026-08-12:

- The WSL distribution is **Ubuntu-22.04**, not Ubuntu 24.04. This matters
  operationally: `watch.ps1` invokes `wsl -d Ubuntu`, which returns
  `WSL_E_DISTRO_NOT_FOUND` against the installed distro name.
- Free space on `D:` is **336 GB of 932 GB**, not 954 GB. Still ample, but see
  Sec. 10.5 — the v2 read batch needs ~50 GB of WSL scratch.

Host: 24 cores, 6.6 GB RAM by default (raised to 9 GB for WSL, see Sec. 8).
Badread is single-threaded (parallelise across jobs); Flye is memory-bound
(limit concurrency).

---

## 3. Plasmid reference set fixed **[DONE]**

`01_filter_plasmids.py` streams the 6.92 GB raw file once and writes a manifest.
Filters fixed before any selection.

| outcome | records |
|---|---|
| scanned | 72,556 |
| **kept** | **61,285** |
| rejected: not "complete sequence" | 6,057 |
| rejected: duplicate sequence | 3,530 |
| rejected: **≥500 kb (chromosome-sized)** | **1,505** |
| rejected: <1 kb | 89 |
| rejected: >0.1 % ambiguous bases | 88 |
| rejected: not declared plasmid | 2 |

The ≥500 kb cut is exactly the proposal's own appendix rule; without it, Badread
would have simulated chromosomes as plasmids.

**Correction to an earlier assumption:** the worry that only 4 of the 10 E. coli
genomes are usable does not apply. Because isolates are *constructed*, a genome's
own plasmid complement is discarded either way — only the chromosome is taken,
and plasmids come from the filtered pool. **All 10 are usable backbones**
(4.76–5.31 Mb).

---

## 4. Simulated cohort v1 **[SUPERSEDED — see Sec. 10]**

`02_build_sim_refs.py`, seed 20260811. Recorded because the v2 rebuild is only
intelligible against what it replaces.

- Plasmid pool: Enterobacterales hosts only; candidates small 7,705 / medium
  11,675 / large 10,906; **240 drawn** (80 per stratum) and extracted.
- **12 isolates, 42 plasmid label units** — strata: 21 small / 13 medium / 8 large.
- Copy number assigned by stratum (small 10–40, medium 2–6, large 1–2).
- Chemistry arms: 6 × `r10_hq` (Badread `nanopore2023` → Flye `--nano-hq`),
  6 × `r9_raw` (`nanopore2020` → `--nano-raw`). Models pinned per arm, never
  defaulted — the default would have made both arms R10.4.1-like and erased the
  chemistry contrast.

### 4.1 The design problem found, and the fix that was applied

The first assembly (sim01 at 50×, default 15 kb fragments) recovered the
chromosome and the 85 kb plasmid but **dropped both small plasmids**. If small
plasmids always drop, plasmid length predicts the label perfectly, B0 saturates,
and the entire ladder becomes uninformative.

Measured the prep effect directly — realised relative depth for sim01's plasmids,
5× probe, fragment length fixed at 6000,4000:

| plasmid | nominal copy no. | rapid (no bias) | ligation (`--small_plasmid_bias`) | effect |
|---|---|---|---|---|
| 2.76 kb | 32 | 27.2× | **1.9×** | **17× depletion** |
| 7.68 kb | 31 | 27.8× | 13.9× | ~2× |
| 84.9 kb | 3 | 3.0× | 2.6× | none |

That is Wick's size-dependent ligation signature reproduced at **zero
bandwidth**, and it is a genuine capability the proposal does not currently
claim (Table I assigns the prep contrast solely to the Wick arm, whose 8–30 GB
transfer is the blocking unknown).

**What this fix got wrong**, discovered only after the full cohort ran: the 17×
depletion is real, but it is measured at the **depth** level and it never reached
the **label** level. At 6000,4000 the plasmids it acts on are shorter than nearly
every read, so Flye reads the circular wrap-around as a repeat and drops them in
*all four cells* regardless of prep or depth. The intervention moved a quantity
that the outcome was already saturated against. See Sec. 10.2.

### 4.2 v1 cohort design

**12 isolates × {ligation, rapid} × {30×, 15×} = 48 assemblies, 168 plasmid
label units.** Fragment length held **fixed** across both preps so prep is not
confounded with N50. Depth is exact, set directly, not subsampled.

> **Correction to how this was reported.** "168 label units" is not a sample
> size. Leave-isolate-out CV and the paired bootstrap both resample **isolates**;
> the 168 units are 42 plasmids × 4 prep×depth cells nested inside 12 isolates
> and carry no independent degrees of freedom. **n = 12.** Quoting 168 overstated
> the evidence by 14×. `06_baselines.py` now leads every printed count with
> isolates and says so explicitly.

---

## 5. Labeller **[DONE, unit-tested]**

`04_label.py` implements Sec. III literally: plasmid is query, GFA segments are
target, `minimap2 -x asm5`, secondary alignments retained, identity applied
**per alignment block before merging**, intervals projected onto the plasmid and
**merged** (never summed).

**Unit test (`test_labeller.py`) — 4 synthetic cases with known answers, all PASS:**

| case | expected | got |
|---|---|---|
| own-size segment | recovered | recovered |
| split across two segments | fragmented | fragmented |
| 4 kb plasmid inside 30 kb segment | absorbed | absorbed |
| not in graph | absent | absent |

Graph-head count and node-positive assignment also correct; the self-loop in the
test GFA is correctly dropped.

**Ambiguity found in the proposal:** precedence between *recovered* and
*absorbed* is undefined — a plasmid can satisfy both (≥95 % into one contig that
is also ≥3× its length). Implemented as absorbed-refines-recovered and flagged
for confirmation.

This file is unchanged by the v2 rebuild. It is validated and correct.

---

## 6. Features and PyG export **[DONE]**

`05_features.py` implements the six node features plus the appendix rules:
chromosome set = segments ≥500 kb with **length-weighted median** depth and
**length-weighted mean** GC; fallback flagged when no segment reaches 500 kb;
degree on the undirected **simple** graph with **self-loops dropped**; self-loop
count kept as a manifest diagnostic and never as a feature.

`05b_export_pyg.py` builds PyG `Data` with the **virtual global node** appended
as the last index and flagged by an extra feature column. S(p) is stored as a
pooling selector, never as a node feature.

Verified on the graphs available:

| tag | segments | simple edges | self-loops | components | isolated | chrom set |
|---|---|---|---|---|---|---|
| sim01_d50 | 2 | 0 | 4 | 2 | 2 | 1 |
| sim01_lig_d15 | 17 | 6 | 4 | 11 | 10 | 4 |

PyG: sim01_lig_d15 → 18 nodes (17 + global), 46 edges (6 simple × 2 + 17 × 2). ✅

Both files are unchanged by the v2 rebuild. **But see Sec. 12.5** — the way the
pooled S(p) vector is *consumed* by the baseline ladder turns out to restate the
label rules, which is a formulation problem sitting downstream of these two
files rather than inside them.

---

## 7. First real results (pilot observations) **[DONE]**

### sim01 at 50× (ligation-like, default 15 kb fragments)

Chromosome recovered (4,758,623 bp, cov 47, circular). The 85 kb plasmid
recovered at **cov 165 vs chromosome 47 = 3.5×** against a nominal copy number
of 3 — confirming Badread's `depth=` header is a per-base coverage multiplier.
Both small plasmids absent. Labels: `{recovered: 1, absent: 2}`.

### sim01 at 15× (ligation) — all three mechanisms in one graph

```
contig_1   7,684 bp  cov 149  circular  -> plasmid_2  RECOVERED
contig_17 48,503 bp  cov  43            \
contig_3  36,394 bp  cov  47            /  sum 84,897 ~= 84,890 -> plasmid_3  FRAGMENTED
(no contig)                                            plasmid_1 (2,760 bp)  ABSENT
chromosome fragmented into 4 contigs >=500 kb (1.98 Mb, 1.06 Mb, 1.06 Mb, 639 kb)
```

Labeller output: `{recovered: 1, fragmented: 1, absent: 1}` — matches manual
reading exactly. The chromosome reference set still resolves without triggering
the fallback rule. Flye at 15× took **96 s**.

This is the mechanistic decomposition the proposal argues for, appearing
naturally in a single simulated isolate — but it did **not** generalise across
the cohort (Sec. 9.1), which is exactly what the repeat-structure work in
Sec. 10.3 is for.

---

## 8. Infrastructure constraint found: the R9.4.1 arm needs ~9 GB **[RESOLVED]**

The assembly batch died three times, each time taking the whole WSL VM with it
and surfacing as an external `SIGTERM` (exit 15). Worth recording because the
diagnosis was initially wrong and because the constraint carries over to the
real cohorts.

**First hypothesis (wrong):** `micromamba run` forwards signals to its process
group, and the manifest line `$(micromamba run flye --version | head -1)` closes
the pipe early. The scripts were rewritten to activate the env by `PATH`
instead. Worth keeping — it is faster and removes a real signal hazard — but it
was **not** the cause: the next run died the same way.

**What actually happened.** A memory sampler at 10 s intervals reported a
comfortable 5.6 GB free throughout, because the spike happened *between
samples*. Re-sampling at **1 s** showed available memory draining 6,212 MB → 14 MB
and the kernel OOM-killing the VM. The "external SIGTERM" was the VM dying.

**The discriminating test:**

| isolate | arm | Flye mode | result | min free RAM |
|---|---|---|---|---|
| sim03 | r10_hq | `--nano-hq` | ok, 71 s, 40 segments | 4,701 MB |
| sim02 | r9_raw | `--nano-raw` | OOM, VM killed | **0 MB** |

Badread's `nanopore2020` model yields ~91 % identity reads (observed Flye
alignment error rate 0.0896), and Flye's k-mer index on them is far larger than
on `nanopore2023`. `--asm-coverage 20 --genome-size 5m` did **not** prevent the
OOM.

**Resolution.** WSL's default allocation is 50 % of host (6.8 GB of 13.7 GB).
`~/.wslconfig` now sets `memory=9GB`, `swap=12GB`; sim02 then assembled in 120 s
— but with only **4 MB** to spare, so swap alone is not a safe margin. Since
Flye's peak scales with thread count, threads are pinned per arm:
**6 for `r9_raw`, 10 for `r10_hq`**.

**Carry-over to the real cohorts — this is not a local quirk.** Wick is R9.4.1
and NEKSUS is R10.4.1, so roughly half the eventual corpus assembles under the
expensive regime. The Sec. IV compute estimate ("trains in minutes… memory well
under 1 GB") is about the *GNN* and remains true; the *assembly* side needs
≥9 GB per concurrent Flye job on R9-era reads. Any plan that assumes assemblies
can be parallelised across cores will be bound by memory, not cores: on a
16 GB machine that is one job at a time, not sixteen.

**Re-verified 2026-08-12 for the v2 batch.** `08_batch_flye.sh` is still strictly
sequential — no `&`, no job pool — and thread pinning is intact at 6/10. A dead
global `THREADS` variable was removed, and `10_run_rest.sh` no longer passes
`JOBS=2 THREADS=11` to it: both were inert (`08` reads neither) but they read as
an instruction to parallelise the arm that OOM-killed the VM three times.

---

## 8b. Assembly batch v1 **[DONE — corpus now superseded]**

**49 of 49 graphs built, zero Flye failures, zero read sets left.**
48 cohort cells + the `sim01_d50` pilot. Reads deleted after graph extraction.

Flye cost by arm (mean): `--nano-hq` ≈ 80 s, `--nano-raw` ≈ 190 s. The R9 arm
runs 2–3× slower as well as needing ~2× the memory.

---

## 9. RESULTS from cohort v1 **[SUPERSEDED — do not quote]**

These numbers validated the *pipeline*. They do not characterise the hypothesis,
and they will be discarded at the v2 rebuild. Retained here because the v2 design
decisions are only justifiable against them.

### 9.1 Realised class distribution

All 49 graphs labelled; every designed plasmid labelled exactly once, none
invented, none lost.

| class | n | % |
|---|---|---|
| recovered | 119 | 69.6 % |
| absent | 32 | 18.7 % |
| fragmented | 8 | 4.7 % |
| absorbed | 1 | 0.6 % |
| grey 50–95 % band (held out) | 11 | 6.4 % |
| **total** | **171** | |

**The four-class decomposition did not materialise.** Fragmented 8 and absorbed
1 are the finding. A mechanism observed once has not been tested; it has been
glimpsed. Sec. 10.3 explains why — and it is a property of the cohort
construction, not of the labeller.

Failure rate was reported as **26.2 / 21.4 / 21.4 / 23.8 %** across (lig,15×),
(lig,30×), (rap,15×), (rap,30×) — "broadly flat".

> **Recomputed 2026-08-12 under the v2 stratum definitions**, over
> *interventional* units only (i.e. excluding the new <4 kb band):
> **9.7 / 3.2 / 6.5 / 9.7 %**, a range of 6.5 percentage points.
>
> The apparent ~23 % failure rate was almost entirely sub-4 kb plasmids:
> **31 of the 32 `absent` units sit below 4 kb.** Of 126 interventional units
> there is exactly **one** absent. The positive class was a constant wearing an
> observation's clothes, and that single fact explains the flat cells, the
> absent-dominated positive class, and the `is_empty` degeneracy in 9.5 below.

### 9.2 Labels track the physical cause

Median realised relative depth by assigned label: **absent 2.71 · fragmented
3.01 · recovered 4.90 · absorbed 32.40**. Absorbed sitting at 32× is the
expected signature — a high-copy plasmid swallowed into a much longer contig.

Failure rate by size band (old bounds): **small 44.0 % · medium 23.1 % · large
3.1 %** — a monotone gradient. Under the v2 bounds the 44 % small figure splits
into a near-deterministic tiny band and a much milder 4–12 kb band.

**Prep and depth each flip 5 plasmid-pairs.** Under the stricter accounting in
`06_baselines.py` (interventional units, prep-or-depth, class-level) the count is
**7 discordant plasmids over 31 units, of which 4 are recovered↔failed flips**.
Against a floor of 20, that is why the matched-pair test could not run.

### 9.3 Graph statistics

49 graphs · segments **2–38** (median 18) · simple edges **0–45** (median 8) ·
isolated nodes median 7, max 18 · **326 self-loops** total (dropped from degree,
kept as diagnostic) · **zero** isolates hit the fragmented-chromosome fallback.

Simulated vs Wick reference plasmid sizes: median **9,546 bp** vs **10,719 bp**;
<10 kb fraction **21/42** vs **14/29**. Reassuringly close, though this compares
*reference* size distributions, not graph statistics — the real
simulated-vs-real graph comparison still needs real assemblies.

### 9.4 The baseline ladder

Leave-isolate-out CV grouped by isolate, PR-AUC, paired bootstrap over isolates.
**n = 12 isolates.** Prevalence 0.248 over 157 nested units.

| model | Arm A (`is_empty` in) | Arm B (`is_empty` out) | Arm C (frag/absorbed only) |
|---|---|---|---|
| B0 "length only" | 0.853 | 0.867 | 0.805 |
| **B0+ primary null** | **0.930** | **0.932** | **0.367** |
| B1 + hand-engineered topology | 0.975 | 0.977 | 0.741 |
| B2 gradient boosting | 0.942 | 0.942 | 0.710 |
| *`is_empty` ALONE* | *0.826* | — | — |

Paired bootstrap, 95 % CI on the PR-AUC difference:

| comparison | Arm A | Arm B | v2 verdict |
|---|---|---|---|
| B1 − B0+ | +0.046 [+0.002, +0.136] | +0.047 [−0.000, +0.142] | **UNRESOLVED** |
| B2 − B0+ | +0.008 [−0.041, +0.082] | +0.008 [−0.052, +0.091] | fails in both arms |
| B0+ − B0 | +0.073 [−0.011, +0.175] | +0.063 [−0.017, +0.161] | fails in both arms |

B1 − B0+ clears in one arm and not the other, on an axis (`is_empty` in vs out)
that should not decide the answer. `06_baselines.py` now labels exactly this
pattern **UNRESOLVED** rather than letting it be reported as clearing or failing.

> **Correction: "B0 length only" is not what was fitted.** `cols_for(["log10_length"])`
> selects columns 0 and 7 of the pooled vector — the **mean and max log length of
> the supporting contigs** — not the plasmid's own length. The plasmid length is
> loaded into the row dict and never enters the design matrix. Measured on the
> same data: the model as coded scores **0.867**; a genuine plasmid-length-only
> baseline scores **0.815**. This is a formulation issue, described in full at
> Sec. 12.5, and it is *not* fixed by the v2 cohort rebuild.

### 9.5 What these numbers do and do not say

**The `is_empty` degeneracy is real and large.** A single binary indicator
reaches **PR-AUC 0.826** against a prevalence of 0.248. Sec. III defines
S(p) = ∅ ⟺ absent, so this is definitional, not learned. Arm B confirms that
deleting the *column* does not fix it: scores are unchanged (0.977 vs 0.975),
because an absent plasmid still carries an all-zero pooled vector. Verified
exactly on this data: **30 of 30 absent units have an entirely zero pooled
vector, and no non-absent unit does.**

**Arm C is underpowered and says nothing.** With absence removed there are only
**9 positives**; B0+ scores *below* B0 and the intervals span ±0.8.

**None of this tests the hypothesis.** The claim is about *learned relational
topology*, which is B2 → B4. B1 is hand-engineered local topology. B3/B4 do not
exist yet. Further, the positive rate (24.8 %) is roughly double the real 9–12 %
because the cohort is 50 % small plasmids against Johnson's 33 % and the absent
class is dominated by plasmids below ~4 kb. **These are pipeline-validation
numbers, not scientific findings.**

---

## 10. Cohort v2 — the rebuild **[RAN — see §16; superseded by v3, §17]**

`02_build_sim_refs.py` rewritten. The pipeline around it is untouched:
`01_filter_plasmids.py`, `04_label.py`, `05_features.py`, `05b_export_pyg.py`
and `test_labeller.py` are all unchanged, because none of them was at fault.

### 10.1 Scale — 12 → 50 isolates

**50 isolates × {ligation, rapid} × {30×, 15×} = 200 assemblies.** 25
compositions cycled twice give **174 plasmids / 696 label units**, of which ~640
are interventional.

The reason is not "more data is better". It is that n was 12 for every test that
mattered, and B1 − B0+ landed at +0.002 in one arm and −0.000 in the other. That
is the signature of the underpowered regime, not of a marginal effect.

Pool raised to 200 accessions per stratum (800 total). Candidates available after
filtering: tiny 3,175 / small 4,990 / medium 11,215 / large 10,906, against a
draw of 28 / 64 / 54 / 28 — 3–7× headroom, which the draw needs because it
filters candidates on whether they stay in band after repeat injection.

### 10.2 Strata — recentred, plus a separate deterministic band

| stratum | bounds | copy number | analysis set |
|---|---|---|---|
| **tiny** | 1,000 – 4,000 | 10–40 | `deterministic_tiny` — **excluded**, reported apart |
| small | 4,000 – 12,000 | 8–30 | interventional |
| medium | 12,000 – 100,000 | 2–6 | interventional |
| large | 100,000 – 500,000 | 1–2 | interventional |

Sub-4 kb is a *different* failure mode, not a more severe version of the same
one: the plasmid is shorter than nearly every read, the wrap-around reads as a
repeat, and Flye drops it in all four cells regardless of prep or depth. Keeping
it, labelling it and excluding it is more honest than either deleting it or
letting it dominate the positive class — which is what it did in v1.

**Fragment length: 6000,4000 → 15000,13000**, identical across both preps so
prep never confounds with N50. Depletion operates around the fragment length, so
raising it moves the affected band up onto 4–12 kb, where dropout is
probabilistic and where the reported ligation depletion actually operates.

> **Stated concern.** Raising the mean also raises wrap-around for the 4–12 kb
> band itself — the same mechanism that made sub-4 kb deterministic at the old
> setting. The wide sd (13000) is what should keep it probabilistic: fragments
> span roughly 2–28 kb, so a 4–12 kb plasmid sees a mix rather than uniformly
> longer reads. This is **not assumed**. If the band goes deterministic again,
> the cell-spread gate and the discordance gate both fire before anything is
> scored.

### 10.3 Repeat structure — the new capability

v1 drew plasmids independently from a pool, so they shared no sequence with each
other or with the chromosome. Fragmentation is repeat-mediated and absorption
needs homology to a longer contig, so **nothing in the cohort could drive either
mechanism**. That is why fragmented came out 8 and absorbed 1.

This is the one axis where *constructed* isolates give leverage a real cohort
cannot: the mechanisms can be induced on demand and audited exactly.

**20 of 50 isolates** (40 %, balanced 10/10 across the two chemistry arms) carry
one IS-like family: a 1–2 kb segment lifted from that isolate's own chromosome —
so its base composition is genuinely genomic and it has a real native locus —
re-inserted at ~99 % identity, substitutions only, in one of two roles:

| role | n | placement | targets |
|---|---|---|---|
| `frag` | 10 | 2–4 copies in the chromosome, 2–3 in one mid-size plasmid, 1–2 in a second plasmid of the same isolate | 20 carrier plasmids → **80 label units** |
| `absorb` | 10 | IS-mediated cointegrate: a ~99 % copy of the whole plasmid embedded in the chromosome between two IS copies, free plasmid dropped to copy number 1 | 10 plasmids → **40 label units** |

The floors are 20 each, so these need conversion rates of 25 % and 50 %
respectively. `09_report.py` §1b prints both conversion rates, so a mechanism
that under-fires is visible immediately and the knob to turn is identified.

> **Deviation from the brief, stated not buried.** The specification asked for
> 1–2 kb shared segments only. Dispersed 1–2 kb elements **cannot** produce
> `absorbed`: that rule needs ≥95 % of the *plasmid* inside one contig, and only
> the shared bases align — a 2 kb element in a 5 kb plasmid is 40 % coverage, not
> 95 %. Reaching 95 % would require a plasmid that is a tandem array of ~20
> copies, which is not IS-like biology. The cointegrate is the mechanism that
> actually produces absorption in nature and it *is* IS-mediated, so it is
> implemented as the `absorb` role and every such placement is flagged
> `kind=cointegrate` in `designs.json`.

### 10.4 Provenance — the cohort now replays from one seed

- `designs.json` becomes `{"provenance": {...}, "isolates": [...]}` carrying the
  seed, the isolate-seed stride, all four stratum bounds, copy ranges, every
  repeat parameter, all 25 compositions, and the input paths; then per isolate:
  its own seed, composition, chromosome, plasmid accessions, and **every repeat
  placement with both original and final byte offsets**.
- `ground_truth.csv` gains nine columns: `plasmid_len_original`,
  `stratum_designed`, `analysis_set`, `repeat_family`, `repeat_role`,
  `repeat_copies_in_plasmid`, `repeat_bases_in_plasmid`, `repeat_frac`,
  `cointegrated`, `isolate_seed`.
- `07_batch_reads.sh` derives a distinct, collision-free seed per
  (isolate, prep, depth) cell and writes it plus the error model to scratch;
  `08_batch_flye.sh` carries both into `manifest.json`. Reads are deleted after
  graph extraction, so unrecorded means unreproducible.

Every placement was verified byte-for-byte against the written FASTA. **That
check found a real bug**: in the `absorb` role the chromosome received two
rounds of insertions, and the second silently invalidated the offsets recorded
by the first (verified identity 0.26 against a recorded 0.99). All insertions
into one replicon now happen in a single pass, and the function's docstring says
why so a future edit cannot reintroduce it.

### 10.5 Operational notes for the run

**Pre-flight is destructive and mandatory.** Every isolate's plasmid complement
changes. `04_label.py` and `05_features.py` walk *every* directory under `asm/`
and label it against `sim/<isolate>_ref.fasta`, and nothing downstream is keyed
to a cohort version — so a surviving v1 `asm/sim01_lig_d15` would be labelled
against plasmids it never contained, come back all-`absent`, and flow into the
feature matrix as real data. `02_build_sim_refs.py` now **refuses to start**
while that state exists (`check_stale_cohort()`, override `REBUILD_FORCE=1`) and
prints the commands:

```
rm -rf work/asm work/graphs
rm -f  work/sim/sim*_ref.fasta
rm -f  work/results/{node_features.csv,graph_diagnostics.json,pyg_dataset.pt,label_validation.csv}
bash   work/scripts/reset_scratch.sh
```

**Capacity.** Measured from the v1 manifests: ~1.0 GB of FASTQ per isolate across
its four cells. `07_batch_reads.sh` generates all 200 read sets before
`08_batch_flye.sh` starts, so peak WSL scratch is **~50 GB** (v1 was ~12 GB)
against 336 GB free on `D:`. It fits, but consider slicing Phase A into two
passes of 25 isolates and interleaving Phase B rather than staging all 50.

**Runtime.** 200 assemblies, 100 on the R9 arm at ~190 s each and ~9 GB RAM each,
100 on the R10 arm at ~80 s. Sequential, so roughly 7.5 hours of Flye alone.

---

## 11. The gates **[SUPERSEDED — now THREE gates, all passing; see §16.1 and §17.1]**

> **This section describes the four-gate panel as it stood for v1 and v2.** The
> `fragmented ≥ 20` gate was **withdrawn** on 2026-08-14 after three probes
> showed the class is not inducible in this simulator (0/164) — see §16.1. The
> remaining three all pass on v3 (§17.1). The text below is kept as the record
> of what was required at the time, not as the current criteria.

`09_report.py` §0b prints this panel; `06_baselines.py` enforces the first one
and **exits without analysing** if it fails, writing `baselines_gate_failed.json`
rather than clobbering the previous `baselines.json`.

| gate | floor | cohort v1 | status |
|---|---|---|---|
| discordant plasmids (prep or depth) | ≥ 20 | 7 (4 binary flips / 31 units) | **FAIL** |
| fragmented units, interventional | ≥ 20 | 8 | **FAIL** |
| absorbed units, interventional | ≥ 20 | 1 | **FAIL** |
| cell failure rates spread | ≥ 10 pp range | 6.5 pp | **FAIL** |

All four fail on v1. That is the finding, not a fault.

The discordance gate matters most because the matched-pair comparison is this
arm's primary contribution: same plasmid, same isolate, **identical length**, one
prep drops it and the other does not. Inside a matched pair B0 is at chance *by
construction*, so anything beating chance is provably using something beyond
length. At 7 pairs that test cannot run.

**B3/B4 (the GNN) are not run until all four pass.** Numbers from a cohort that
fails these gates are discarded at the next rebuild anyway.

---

## 12. Audit of the working directory, 2026-08-12 **[DONE]**

### 12.1 Divergences between this document and the code

Six, all now corrected in place above: the B0 mislabel (§9.4), WSL Ubuntu 24.04
vs 22.04 (§2), 954 GB vs 336 GB free (§2), "168 label units" as a sample size
(§4.2), the flat cell-failure gradient being a tiny-plasmid artefact (§9.1), and
an appendix file map missing `11_assemble_all.sh` and the twelve helper scripts.

### 12.2 What would have broken at 50 isolates

- **Stale `asm/` and `graphs/`** — silent garbage, not an error. Fixed by
  `check_stale_cohort()`. This was the one that would have cost a full run.
- **`04b_validate_labels.py`** hardcoded `("small","medium","large")` in three
  places and would have silently dropped the entire tiny band from all three of
  its reports. Fixed with a `STRATA` constant; the determinism check now runs
  explicitly on `small` (4–12 kb) and says so.
- **`10_run_rest.sh`** passed `JOBS=2 THREADS=11`. Inert, but a parallelisation
  trap. Removed.
- **`watch.ps1`** hardcodes `$total = 49`, calls `wsl -d Ubuntu` (returns
  `WSL_E_DISTRO_NOT_FOUND`), and reads a task-output path under a stale session
  UUID and a different Windows user. Not fixed — monitoring only.
- **`03_sim_assemble.sh`** writes `simNN_dDD` tags that bypass the 2×2 filters in
  `06` and `04b` but still enter `node_features.csv` and the §2 graph
  statistics. That is exactly how the `sim01_d50` pilot got into the v1 corpus.

### 12.3 Old constants still in the tree

`test_prep_arms.sh:28` still probes at `--length 6000,4000` — harmless as a
one-off measurement script, but it no longer matches the cohort and will mislead
anyone re-measuring the prep effect. `03_sim_assemble.sh` passes no `--length`
at all and inherits Badread's default, which happens to *be* 15000,13000, so it
now agrees with the cohort by accident rather than by intent. No old stratum
bounds survive; the `< 10000` tests in `09_report.py` are Johnson's and Wick's
reporting convention for "small plasmid", kept so that comparison stays
commensurable with the published tables and now commented as deliberately *not*
the stratum bound.

### 12.4 The NEKSUS leak guard

`scripts/leak_guard.py` is new. `assert_no_leak_columns()` hard-fails — never
warns — on any column name containing `circular` or `completeness`, in any casing
or separator style. `load_neksus_contigs_summary()` splits those columns off **at
the parser**: there is no code path that returns them inside `features`. Run
against the real 3,917-row file it strips five columns, three more than the two
originally identified. Both are retained as diagnostics. `06_baselines.py`
imports it and asserts the feature list at the point of definition.

### 12.5 Remaining label-leak paths **[OPEN — formulation, not cohort]**

Three beyond the two already known (Flye's circularity flag; `is_empty` / the
zero pooled vector). All three are properties of how the plasmid head consumes
S(p), so **none is fixed by the cohort rebuild** and none has been changed
unilaterally. They are written up as **`PROPOSAL_FIXES.md` Section 1, open item B**
— the one decision still outstanding on the confirmatory axis.

1. **The pooled support-segment length is close to a restatement of the label
   rule.** max(log10 support length) − log10(plasmid length), by class:
   absorbed **+0.488**, recovered **−0.000**, absent **−3.440**. The `absorbed`
   rule's own threshold is log10(3) = **+0.477**. The feature is the decision
   boundary.
2. **B0 carries the absence indicator.** All 30 absent units have an entirely
   zero pooled vector and B0's two columns are part of it, so the supposed
   no-topology floor is contaminated at exactly the place the confirmatory
   comparison runs.
3. **|S(p)| ≥ 2 *is* the fragmented rule**, and mean ≠ max on any pooled feature
   reveals it.

Housekeeping rather than leak: `node_features.csv` ships `_raw_len`,
`_raw_depth`, `_gc`; `plasmid_labels.csv` ships `best_cov_frac`,
`cumulative_cov_frac`, `best_segment_len`. Everything reads them by column name
today, but one change to a `read_csv().values` idiom pulls them straight in.

### 12.6 Provenance gaps that remain

The v2 cohort *does* replay exactly from `SEED`. Three gaps remain outside it:
the cohort manifests record `flye_version` but no Badread or minimap2 version
(`03_sim_assemble.sh` recorded Badread's, but `03` is not the cohort path);
`01_filter_plasmids.py` prints its filter counts to stdout only, persisting no
checksum or record count of the 6.92 GB input; and **the working directory is not
a git repository**, so no output is pinned to the script version that produced it.

---

## 13. Still pending

*Updated 2026-08-15 after the v3 run. Closed items struck through; see §17.*

| item | blocker |
|---|---|
| ~~Fix B0~~ | **DONE** — PROPOSAL_FIXES item A closed; the contamination it removed was +0.31 to +0.67 PR-AUC (§17.3) |
| ~~Run the rebuild~~ | **DONE** — v3, 200/200, zero failures (§17) |
| ~~Score the cohort~~ | **DONE** — all three gates pass, ladder in §17.4 |
| **Real cohorts (Wick first)** | **now the critical path.** §17.6: 79 edges over 200 graphs. B3/B4 cannot be tested on simulated data |
| B3 / B4 (GNN) | **real assembly graphs**, not the gates — those now pass. Building on this corpus would be an MLP |
| Confirmatory-axis decision (open item B) | a decision, not code; needed to *interpret* B4 |
| Correctness gate against Johnson | needs Johnson ONT reads (not on disk) |
| Simulated-vs-real *graph* statistics | needs any real assembly |
| Johnson reference FASTAs via Unicycler | needs Johnson Illumina reads |
| Cointegrate conversion back to ~90 % | optional; §17.7 says how |

---

## 14. Proposal edits required

**`PROPOSAL_FIXES.md` is now a checklist, not an archive.** It has three
sections: **open items** (2), **closed items** (13, as a verification table), and
the **evidence archive** carrying the full text of everything closed. The rule at
the top of that document is that a fix *moves* an item into the closed table — it
never adds a new item — so the open count falls as work is done.

**Two items are open, and neither is fixed by this cohort rebuild.**

- **A — B0 is not a length-only baseline** `[CODE, BLOCKING]`. Affects
  `06_baselines.py`; must land **before the v2 cohort is scored**, because every
  rung is a difference from a floor that is currently contaminated. See §9.4 and
  §12.5 above.
- **B — the confirmatory axis runs on a head that restates its own labels**
  `[DECISION]`. Needs no code until B3/B4 — the recommended resolution moves the
  confirmatory comparison to the node and graph heads, which *are* the GNN. It
  does not block the v2 run; it blocks interpreting it.

The two worth singling out:

**Sec. III states S(p) = ∅ ⟹ the plasmid is absent**, so the `is_empty` indicator
predicts a positive subclass *by definition rather than by learning* — on the
confirmatory axis (B0+ vs B4). Measured: `is_empty` alone reaches PR-AUC 0.826
against a prevalence of 0.248.

**The same class of problem is larger than one indicator.** The pooled S(p)
vector also encodes the `absorbed` threshold and the `fragmented` support-count
rule, and the supposedly intrinsic B0 floor is contaminated by absence. The
plasmid head reads several of its own label rules back out of its features.
This is the formulation problem the cohort rebuild does not touch.

---

## 15. Defensible mid-review claim

*Rewritten 2026-08-15 for the v3 result.*

> The end-to-end read → assembly → labelling → feature → graph pipeline is
> implemented and validated on controlled simulated isolates with known ground
> truth: **449 assemblies across three cohorts with zero assembler failures**, a
> unit-tested labeller, and complete label accounting on every run. Two cohorts
> were rejected by adequacy criteria fixed in advance before any model was
> scored. The third passes all of them — 69 discordant matched pairs, a 13.9 pp
> cell-failure spread against a 10 pp floor, and a realised spread within
> 1.1 pp of the value predicted from a design probe *before* the cohort was
> built. A contaminated baseline was found and corrected: the nominal
> length-only floor was reading absence out of its own pooled feature vector,
> worth **+0.31 to +0.67 PR-AUC**, which had made every previously reported rung
> uninterpretable. Three purpose-built probes then established that
> repeat-mediated fragmentation is **not inducible** in a Badread → Flye
> simulation (0/164 across repeat length, copy count, identity, chemistry and
> depth), so that class was withdrawn from the gate panel on evidence rather
> than relaxed after a near miss. The same probes establish the arm's boundary:
> **79 assembly-graph edges across 200 graphs, 73 % of graphs edgeless**, so the
> relational-topology hypothesis cannot be tested on simulated data and requires
> real assembly graphs.

The claim is **not** "the GNN works" and **not** "topology beats the null" —
B1 − B0+ is formally UNRESOLVED and is reported as such. It is that the
pipeline is validated, two cohorts were rejected by pre-registered criteria,
the third passes them against a stated prediction, a real baseline leak was
found and quantified, and the arm's limits were established by experiment
rather than assumed.

---

## 16. Cohort v2 ran, failed 2 of 4 gates, and was diagnosed **[DONE, 2026-08-14]**

**200/200 assemblies, zero Flye failures, 11.6 h.** Archived to `cohort_v2/`
with its own README. Gate panel: discordant **27** ✅ · absorbed **40** ✅ ·
fragmented **19** ❌ · cell spread **5.5 pp** ❌.

### 16.1 Three probes, run before touching the builder

Rather than rewrite the cohort design on a theory and spend another 12 h, each
question was settled on one purpose-built isolate. All three live in separate
roots (`probe/`, `probe2/`, `probe3/`) and never touched `work/asm`.

**Probe 1 — where does prep change the label?** Eight plasmids at 5–12 kb, copy
number pinned at 15× so length is the only variable. Ligation depletes ~10× at
the *depth* level across the whole ladder, but the *label* only flips between
**9,083 and 10,046 bp**. Below that both preps fail; above it both succeed.

**Probes 1+2 — can fragmentation be induced?** No.

| axis | tested | fragmented |
|---|---|---|
| repeat length | 2, 3, 5, 8, 15, 25 kb | 0 |
| copies in plasmid / chromosome | 2–6 / 0–4 | 0 |
| identity | 0.990, 1.000 (verified verbatim) | 0 |
| chemistry × depth | r10_hq/r9_raw × 30×/15× | 0 |

**0/164 including the cohort's own 0/148.** Flye assembled every carrier
complete and circular and emitted the repeats as separate high-coverage contigs
beside them — it detects them and resolves them anyway. Depth is not the lever
either: v2's medium/large units below 1.5× relative depth came out **56
recovered vs 6 fragmented**, and fragmented sits at a *higher* median relative
depth than recovered (3.99 vs 2.87).

**Probe 3 — length × copy number over all four gate cells.** 27 plasmids,
7.6–11.9 kb, crossed against copy {8, 15, 25}, plus an r9_raw check.

| band | n | L30 | R30 | L15 | R15 | spread |
|---|---|---|---|---|---|---|
| 1,000–8,400 | 5 | 80.0 | 60.0 | 60.0 | 60.0 | 20.0 pp |
| 7,500–12,500 | 27 | 33.3 | 18.5 | 33.3 | 14.8 | 18.5 pp |
| **8,400–10,900** | **19** | **26.3** | **10.5** | **31.6** | **5.3** | **26.3 pp** |

`probe_v3c_report.py` also proposes 8,828–9,251 at a "100 % flip rate"; that is
min/max over the 3 plasmids flipping at both depths and is **overfitted — not
used**. The band chosen rests on 19 plasmids.

### 16.2 The finding that outlives the cohort **[OPEN — scope]**

**45 real edges across 200 graphs · 175 of 200 edgeless · 698 of 768 nodes
isolated (91 %).** The three probes, built specifically to create branch points,
produced **one edge between them**. Badread samples uniformly from a perfect
circular template and Flye traverses that circle correctly whatever is inside it.

**The simulation arm does not generate assembly-graph topology, so it cannot
test the relational-topology hypothesis.** B3/B4 on this corpus would be an MLP
with global pooling. That claim needs real assembly graphs — Wick's 7 isolates
remain the cheapest route.

### 16.3 Cohort v3 **[BUILT, NOT RUN]**

| | v2 | v3 |
|---|---|---|
| `small` band | 4,000–12,000 | **8,400–10,900** (measured) |
| `small` copy | 8–30 | 8–25 (inside what was tested) |
| `tiny` band | 1,000–4,000 | 1,000–8,400 — now *the deterministic band*, key name kept to avoid a 4-file rename |
| `small` share of interventional | 44 % | **57 %** |
| repeat roles | frag 10 / absorb 10 | **absorb 20** — `frag` withdrawn |
| cointegrate slot | shortest (always `small`) | **medium/large** — cointegrates are constants and were eating the prep window |
| gates | 4 | **3** — `fragmented` floor withdrawn on probe evidence |

Realised: 170 plasmids, 680 units, 158 interventional, small 90 / medium 46 /
large 22 / tiny 12, 20 cointegrates (16 medium, 4 large, **0 small**), 0 out of
band, 0 FASTA-vs-ground-truth length mismatches.

> **Prediction, stated before the run:** cell spread **≈ 15 pp** (0.57 × 26.3)
> against a 10 pp floor. Absorbed ≈ 72 units at v2's 90 % conversion. If the
> spread comes in under 10 pp, the band is wrong and the fragment length is the
> next thing to move — not the band.

Bounds are valid for `FRAGLEN 15000,13000` **only**; change it and re-run
`probe_v3c.*` before trusting them.

---

## 17. RESULTS — cohort v3 **[DONE, 2026-08-15]**

**200/200 assemblies, zero Flye failures.** Built across several driver
invocations after a mid-run pause; that does not affect the corpus, because
every cell's Badread seed derives from the isolate's index in `designs.json`
rather than from run history. Labelling accounting is clean: **680 label units
over 200 graphs, every designed plasmid labelled exactly once, none invented,
none lost.**

### 17.1 All three gates pass

| gate | floor | v1 | v2 | **v3** |
|---|---|---|---|---|
| discordant plasmids (prep or depth) | ≥ 20 | 7 | 27 | **69** (56 binary flips) |
| absorbed units, interventional | ≥ 20 | 1 | 40 | **47** |
| cell failure-rate spread | ≥ 10 pp | 6.5 | 5.5 | **13.9 pp** |
| ~~fragmented units~~ | — | 8 | 19 | 17 · **withdrawn**, see §16.1 |

**The prediction held.** §16.3 stated ≈ 15 pp before the cohort was built,
from 0.57 × 26.3 pp measured in probe 3. Realised **13.9 pp**.

```
prep=lig  depth=15x   failure 25.3%        prep=rap  depth=15x   failure 13.9%
prep=lig  depth=30x   failure 22.8%        prep=rap  depth=30x   failure 11.4%
```

Prep separates the cells (lig ≈ 24 %, rap ≈ 13 %) with a depth gradient inside
each — the designed structure, not a lucky margin.

### 17.2 Class distribution and the interventional axes

| class | n | % | interventional |
|---|---|---|---|
| recovered | 471 | 69.3 | 471 |
| absent | 98 | 14.4 | 52 |
| absorbed | 47 | 6.9 | 47 |
| gray 50–95 % (held out) | 47 | 6.9 | 45 |
| fragmented | 17 | 2.5 | 17 |
| **total** | **680** | | **632** |

By stratum: tiny 95.8 % failure (48 units, excluded) · **small 17.2 %** ·
medium 26.1 % · large 6.8 %. The determinism check passes on `small`.

**Prep flips 68 plasmid-pairs** (v2: 26), asymmetric in the designed direction —
53 lig-fails/rap-ok against 15 the other way. **Depth rescues 28** (v2: 7).

Dose-response by median realised relative depth: absorbed 0.97 · absent 1.72 ·
gray 2.57 · recovered 2.98 · fragmented 3.95. Absorbed sitting lowest is the
cointegrate signature — the free plasmid is dropped to copy number 1.

### 17.3 The B0 fix was not cosmetic

`B0_pooled − B0` measures the contamination directly, and it is large in every
arm: **+0.306 [+0.231, +0.383]** (A), **+0.334 [+0.259, +0.424]** (B),
**+0.667 [+0.559, +0.771]** (C). All three clear.

In Arm C the true length-only floor scores **0.135 against a prevalence of
0.120** — chance — while the old pooled definition scored **0.806**. Every rung
in every earlier report was measured from a floor sitting 0.3–0.67 PR-AUC too
high, purely on pooled-S(p) contamination. `is_empty` alone is now **0.504**,
down from v1's 0.826.

### 17.4 The ladder

Leave-isolate-out CV, **n = 50 ISOLATES**, PR-AUC, paired bootstrap over
isolates. 587 nested units / 116 positives (A, B); 535 / 64 (C).

| model | Arm A (`is_empty` in) | Arm B (out) | Arm C (frag/absorbed only) |
|---|---|---|---|
| **B0** plasmid's own length | 0.553 | 0.239 | **0.135** |
| B0_pooled *(the old B0)* | 0.864 | 0.573 | 0.806 |
| **B0+** primary null | 0.891 | 0.617 | 0.846 |
| B1 + hand-engineered topology | 0.892 | 0.865 | 0.800 |
| B2 gradient boosting | **0.956** | **0.956** | **0.947** |
| *`is_empty` ALONE* | *0.504* | — | — |

| comparison | Arm A | Arm B | cross-arm verdict |
|---|---|---|---|
| B1 − B0+ | +0.001 [−0.007, +0.009] | +0.245 [+0.165, +0.321] | **UNRESOLVED** (ci_lo changes sign) |
| B2 − B0+ | +0.062 [+0.001, +0.137] | +0.331 [+0.243, +0.417] | clears in both |
| B0+ − B0 | +0.338 [+0.234, +0.435] | +0.379 [+0.312, +0.461] | clears in both |
| B0_pooled − B0 | +0.306 [+0.231, +0.383] | +0.334 [+0.259, +0.424] | clears in both |

**B2 − B0+ clears in both arms** — nonlinearity buys something here, unlike v1
where it bought nothing (+0.008, interval spanning zero).

**B1 − B0+ is UNRESOLVED and must not be reported either way.** +0.001 in one
arm and +0.245 in the other, on an arm that should not decide the answer.

Arm T (<4 kb) was skipped: 46 units, all positive, nothing to fit.

### 17.5 What these numbers do and do not say

**They do say** the cohort is adequate for what this arm claims — the pipeline
is validated end to end on 200 assemblies with zero failures, the prep and
depth contrasts move the label as designed, `absorbed` is inducible, and the
baseline floor is now honest.

**They do not test the hypothesis.** B2 is nonlinearity on hand-engineered
features, not learned relational topology. B3/B4 still do not exist, and §17.6
is why building them on this corpus would be pointless.

### 17.6 Topology: better, still not enough **[OPEN]**

| | v2 | **v3** |
|---|---|---|
| real edges over 200 graphs | 45 | **79** |
| graphs with zero edges | 175 | **146** |
| isolated nodes | 91 % | **86 %** |
| max edges in one graph | 5 | 8 |

Improved, and still far short. 79 edges across 200 graphs, with 73 % of graphs
edgeless, does not support a message-passing model — B3/B4 here would be an MLP
with global pooling. **§16.2 stands: the relational-topology claim needs real
assembly graphs.** Wick's 7 isolates remain the cheapest route.

### 17.7 One cost of a v3 change, recorded

Cointegrate → absorbed conversion fell to **43/80 = 54 %** from v2's 90 %.
Moving the cointegrate off `small` onto `medium` (§16.3) protected the prep
window, but a medium plasmid at copy number 1 still assembles on its own, so 34
of 80 came back `recovered` instead of collapsing into the chromosome. Net
absorbed is still up (47 vs 40) because there are twice as many cointegrates.
If the rate matters more than the count, cointegrate the *smallest* medium
slots rather than any medium slot.

---

## Appendix: file map

```
work/
  scripts/
    01_filter_plasmids.py    6.9 GB raw query -> filtered manifest        [untouched]
    02_build_sim_refs.py     chromosomes + pool + per-isolate refs        [v2 rewrite]
    03_sim_assemble.sh       single isolate: reads -> Flye -> GFA (pilot path)
    04_label.py              GFA + references -> four label classes + S(p) [untouched]
    04b_validate_labels.py   labels vs ground truth (accounting, dose-response,
                             prep effect, length-determinism check)       [strata fix]
    05_features.py           GFA -> six node features + diagnostics        [untouched]
    05b_export_pyg.py        -> PyG Data with virtual global node          [untouched]
    06_baselines.py          B0/B0+/B1/B2, leave-isolate-out, paired bootstrap,
                             adequacy gate, UNRESOLVED verdicts           [v2]
    07_batch_reads.sh        Phase A: parallel read generation             [v2]
    08_batch_flye.sh         Phase B: assembly, mode + threads pinned by arm [v2]
    09_report.py             everything Phase 0 must report + gate panel   [v2]
    10_run_rest.sh           chains Phase A -> Phase B                     [v2]
    11_assemble_all.sh       restartable sequential assembly driver
    leak_guard.py            circular/completeness blocked at the parser   [new]
    test_labeller.py         labeller unit test (4 known-answer cases)     [untouched]
    env.sh                   activate the plasmid env by PATH, not `micromamba run`
    install_tools.sh install_badread.sh envchk.sh sudochk.sh test_env.sh
    status.sh progress.sh diagnose.sh check_reads.sh reset_scratch.sh
    run_label.sh run_test_labeller.sh
    test_flye_lowmem.sh test_flye_persist.sh test_one_flye.sh
    test_prep_arms.sh test_prep_yield.sh badread_depth.sh badread_help.sh
  refs/                      plasmid_manifest.csv, plasmid_pool.fasta
  sim/                       simNN_ref.fasta, ground_truth.csv, designs.json
  asm/<tag>/                 assembly_graph.gfa, assembly_info.txt, manifest.json,
                             replicon_yield.tsv   (reads deleted after extraction)
  graphs/<tag>/              plasmid_labels.csv, node_labels.csv, aln.paf
  results/                   node_features.csv, graph_diagnostics.json,
                             pyg_dataset.pt, label_validation.csv, baselines.json,
                             baselines_gate_failed.json
    probe_v3.py/.sh          probe 1: prep window (ladder A) + repeat length (B)
    probe_v3b.py/.sh         probe 2: repeat identity / copy count / placement
    probe_v3c.py/.sh         probe 3: length x copy number over all 4 gate cells
    probe_v3*_label.sh       label a probe root via PIPE_ROOT
    probe_v3_report.py       ladder A/B readout
    probe_v3c_report.py      length x copy map (its band suggestion is overfitted)
    12_run_cohort.sh         A->B in isolate slices; args: START SLICE PAR
  probe/ probe2/ probe3/     probe roots -- separate from asm/, never labelled
                             into the cohort. See PROGRESS.md sec. 16.1.
  cohort_v2/                 ARCHIVED v2 corpus: 200 assemblies, 3.2 GB, + README
                             explaining which 2 gates it failed and why
  cohort_v1/                 ARCHIVED v1 corpus: asm/ graphs/ sim/ results/ +
                             README.md. 49 assemblies, 786 MB, moved out of
                             asm/ and graphs/ on 2026-08-13 so the v2 rebuild
                             cannot label them against the wrong references.
                             Invisible to the pipeline (discovery is a
                             non-recursive listdir of asm/ and graphs/).
  watch.ps1                  live batch view (stale: hardcodes 49, wrong distro name)
  PROPOSAL_FIXES.md          checklist: 2 open / 13 closed + evidence archive
  PROGRESS.md                this file
```
