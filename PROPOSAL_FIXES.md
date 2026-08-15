# Proposal fixes

**Revision 12 · 2026-08-14 — 1 open, 4 drafted-but-unapplied, 10 closed.**

---

## 1. OPEN

### B. Confirmatory axis runs on a head that restates its own labels · `[DECISION]`

The pooled S(p) vector encodes three of the label rules:

| rule | recoverable statistic | measured |
|---|---|---|
| *absorbed* — contig ≥ 3× plasmid | max(log10 support len) − log10(plasmid len) | **+0.488** vs the rule's own **log10(3) = +0.477** |
| *absent* — S(p) = ∅ | `is_empty`, or equivalently the zero vector | **PR-AUC 0.826** alone, at **prevalence 0.248** |
| *fragmented* — ≥ 2 contigs | mean ≠ max on any pooled feature | **\|S(p)\| ≥ 2 IS the rule** |

The confirmatory comparison (B0+ vs B4) runs on this head and *absent* is one of
its three positive mechanisms, so the headline test can be carried by a subclass
any model gets for free.

**Recommended:** move the confirmatory test to the **node head**
(fragmented/absorbed — every unit has a non-empty support set by construction)
and the **graph head** (absence, by count calibration). Plasmid head stays for
training, reported descriptively with the `is_empty`-alone baseline beside it.
**No code needed until B3/B4 exist.**

**Needs: a decision, recorded here once made.**

---

## 2. PRE-FLIGHT — before running v2

- [x] **A fixed and asserted** in `06_baselines.py` — closed, see §4
- [x] **Destructive cleanup done** (`rm -rf work/asm work/graphs`, etc.) — stale
      v1 assemblies would be labelled against v2 references, return all-absent,
      and enter the feature matrix as plausible data.
      Done 2026-08-13: the v1 corpus was **archived, not deleted**, to
      `work/cohort_v1/` (49 assemblies + graphs + sim + results, 786 MB, with a
      README). It sits outside `asm/` and `graphs/`, both of which are
      discovered by a non-recursive `listdir`, so it is invisible to the
      pipeline and does not trip `check_stale_cohort()`
- [ ] **`03_sim_assemble.sh` outputs tagged `analysis_set=pilot`**, filtered on
      the tag not the name pattern — this is how `sim01_d50` contaminated v1
- [ ] **Locate or create the proposal document**, then apply items 1, 2, 3, 7
      from §4 to it. Drafted text is ready below; application is blocked only on
      the document existing. **This item currently has no owner.**

---

## 3. TO CHECK — only resolvable after the v2 run

Do not mark these done from reading code.

- [ ] **Four gates all pass:** discordant ≥ 20 (v1: **7**), fragmented ≥ 20
      (v1: **8**), absorbed ≥ 20 (v1: **1**), cell spread ≥ 10 pp (v1: **6.5 pp**).
      Read `09_report.py` §0b **first**. If discordance and spread fail together,
      the 4–12 kb band went deterministic — **widen the fragment sd, don't raise
      the mean.**
- [ ] **Mechanism conversion rates hit target:** IS-carrier → fragmented needs
      25 % of 80 candidate units; cointegrate → absorbed needs 50 % of 40.
      Printed in `09_report.py` §1b.
- [ ] **B0 fix left no other model contaminated** by the pooled vector — re-run
      the pooled-vector audit (item 12's method) against **v2** data once
      labelled, not just against v1.
- [x] **Peak scratch stayed inside plan** — resolved by changing the plan, not
      by confirming it. Staging all 200 read sets before Phase B would have
      peaked at ~42 GB (measured: v1 wrote 10 GB for 12 isolates) against
      **40.7 GB free on C:**, where the WSL `ext4.vhdx` lives and which does not
      shrink once grown. `07_batch_reads.sh` gained `ISO_FROM`/`ISO_TO` and
      `12_run_cohort.sh` now interleaves A→B in slices of 10 isolates.
      **Actual peak 8.7 GB**, and the vhdx did not grow at all (15.4 GB before
      and after). Seeds still derive from each isolate's index in the full
      list, so a sliced run reproduces an unsliced one byte-for-byte.

---

## 4. CLOSED / DRAFTED

| # | Issue | Status / closed by | Verify in |
|---|---|---|---|
| 1 | Gate tolerance undefined | **DRAFTED, not applied** | no target document exists |
| 2 | NEKSUS genus split degenerate | **DRAFTED, not applied** | no target document exists |
| 3 | Johnson refs not deposited | **DRAFTED, not applied** | no target document exists |
| 4 | NEKSUS refs are metadata | `leak_guard.py`, 5 cols stripped | `leak_guard.py`; prose addition also drafted, not applied (below) |
| 5 | *absorbed*/*recovered* precedence | absorbed refines recovered | `test_labeller.py` |
| 7 | Simulation transfer wrong ~23× | **DRAFTED, not applied** | no target document exists |
| 8 | Prep contrast assigned Wick only | `--small_plasmid_bias`, **17×** depletion | `02_build_sim_refs.py` |
| 8b | Prep effect must be at label level | FRAGLEN 15000,13000 | `07_batch_reads.sh` |
| 8c | No repeat-mediated mechanism | 20 isolates, IS + cointegrate | `02_build_sim_refs.py` |
| 9 | Assembly assumed memory-cheap | sequential R9 arm, threads 6/10 | `08_batch_flye.sh` |
| 10 | 168 units reported as sample size | isolate-first, UNRESOLVED flag | `06_baselines.py` |
| 11 | <10 kb band mixed 2 failure modes | 4 strata, tiny excluded | `02_build_sim_refs.py` |
| 13 | distribution ≠ adequacy | 4 gates, class floors | `09_report.py`, `06_baselines.py`; prose addition also drafted, not applied (below) |
| A | B0 was not a length-only baseline | `PLASMID_LEN` column outside the pooled block; B0 uses it alone; old rung kept as `B0_pooled`; 4 import-time asserts | `06_baselines.py` |

**14** merged into **8c**, same argument as 8. **6** subsumed by open item **B**.

**A — how it was closed (2026-08-14).** `plasmid_labels.csv` carries
`plasmid_len` as `len(seq)` taken straight from `sim/<isolate>_ref.fasta` by
`04_label.py`, so that is the reference-FASTA source the fix called for; it was
checked against `ground_truth.csv` for all 174 plasmids (0 mismatches, and
`plasmid_len_original` differs for exactly the 20 repeat-carriers, as it should).
The column is appended after `is_empty` at index `PLEN_COL`, outside
`POOLED_COLS`. Four assertions fire at import, before anything is fitted:
B0 selects exactly one column; that column is `PLEN_COL`; it is not inside the
pooled block; and `PLASMID_LEN` is not in `FEATS_ALL`. Both regressions were
tested by patching the bug back in — reverting B0 to `["log10_length"]` and
adding the intrinsic name to `FEATS_ALL` — and both are rejected at import.
Behaviourally verified: an absent unit and a recovered unit of the same length
now carry an identical B0 value, so the floor cannot see the label.

`load()` additionally hard-exits if `plasmid_len` disagrees between
`plasmid_labels.csv` and `ground_truth.csv`, which is what a stale labelling
against a different cohort would look like.

**Note carried forward, not silently absorbed:** B0+ is still the *pooled*
intrinsic triple, so `B0+ − B0` now spans two changes at once — the added depth
and GC features, and the move from an intrinsic column to a pooled vector that
is zero for every absent unit. `B0_pooled − B0` isolates the second, and is
printed with the other comparisons. Whether B0+ should also become intrinsic is
part of open item **B**, not this fix.

---

## Drafted text, pending a destination document

Held verbatim so nothing is lost a second time. Numbers sourced from PROGRESS.md
§1.2 and §3, not reconstructed.

**When applied:** move the row above from *DRAFTED, not applied* to a normal
closed row with the real section reference, and delete the corresponding block
from here.

### Item 1 — target: wherever the correctness-gate tolerance is specified

> The correctness gate is defined against Johnson replicate 1 under the pinned
> mode: **`flye-raw`, replicate 1, target 30 of 33 plasmids recovered, tolerance
> ±1.** The replicate is named because Johnson's own across-replicate range for
> this mode is **29–31** (30/33, 29/33, 31/33) — wider than the tolerance itself,
> so comparing against an unspecified replicate would make the gate
> unfalsifiable. The across-replicate range is reported alongside the gate result
> as context.

### Item 2 — target: wherever the split-by-genus claim is made

> NEKSUS supports exactly one genus contrast: *Escherichia* (**58** isolates)
> against *Klebsiella* (**29**), across the 92 assembled isolates. The remaining
> genera — *Citrobacter* (2), *Enterobacter* (2), *Serratia* (1) — total **five**
> isolates and cannot support a stratified held-out side, so they are assigned to
> the training side and reported separately rather than used for evaluation. The
> genus split is therefore a single-contrast check on one taxonomic boundary, not
> broad unseen-taxa validation, and is described as such wherever it is reported.

### Item 3 — target: wherever Johnson reference provenance is stated

> Johnson reference plasmids are **not deposited as sequence**; the supplement is
> two CSV tables — per-plasmid read statistics and per-plasmid recovery — and
> neither contains sequence. All 14 references are therefore re-derived with a
> pinned Unicycler, which requires the Johnson **Illumina** reads in addition to
> the ONT reads, a cost the ~1.8 GB transfer figure does not cover. The
> re-derivation is added to the Phase 0 critical path. Because every Johnson
> label now depends on our own Unicycler run rather than the authors', the audit
> subset covers Johnson as a matter of course, not as a spot check.

### Item 7 — target: wherever the simulation transfer size is stated

Transfer column: **"~300 MB"** becomes **"~7 GB raw plasmid query, reduced to a
filtered reference set before use"**, plus:

> The plasmid source is filtered before any simulation: records must declare
> `plasmid` and `complete sequence`, fall in **1 kb ≤ L < 500 kb** (matching the
> ≥500 kb chromosome rule), carry **≤0.1 %** ambiguous bases, and be unique by
> accession and by sequence. This removes **11,271 of 72,556** records —
> including **1,505** that are chromosome-sized — leaving 61,285. Filter counts
> are reported in Phase 0.

### Item 4 (prose half) — target: the feature specification section

> Any column expressing an assembler's own verdict on replicon **circularity** or
> **completeness** is excluded from every feature matrix and retained as a
> diagnostic only. This is enforced at the parser rather than by convention: on
> `contigs_summary_sup_cleaned.csv` the guard strips five columns — `circular`,
> `completeness`, `no_circular_chromosomes`, `no_circular_contigs` and
> `chromosome_circularity`.

### Item 13 (prose half) — target: the evaluation / results-reporting section

> The realised class distribution is reported **against explicit floors fixed in
> advance**, and each class is stated as meeting its floor or not: **≥20
> fragmented** and **≥20 absorbed** units, counted over interventional strata
> only. **No per-class performance claim is made for a class below its floor.**
> Two cohort-level criteria are checked before any model is fitted: **≥20
> discordant matched pairs** — plasmids whose label differs across prep or across
> depth at otherwise identical settings — and **≥10 percentage points** of
> variation in failure rate across the prep × depth cells. A cohort failing any
> criterion is rebuilt rather than analysed, and the failure is reported as a
> refusal to analyse, not as a result.
