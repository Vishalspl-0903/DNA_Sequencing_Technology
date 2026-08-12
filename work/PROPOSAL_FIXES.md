# Proposal fixes

**Revision 11 · 2026-08-12 — 2 open, 4 drafted-but-unapplied, 9 closed.**

---

## 1. OPEN

### A. B0 is not a length-only baseline · `[CODE, BLOCKING]`

`cols_for(["log10_length"])` selects the **mean and max log length of the
supporting contigs**, not the plasmid's own length — which never enters the
design matrix at all.

It is also contaminated: all **30 of 30** absent units carry an entirely zero
pooled vector, B0's two columns included, so the model designated as the
no-topology floor already encodes the class it is meant to be ignorant of.
Scored both ways on the same data: **as coded 0.867; true plasmid-length-only
0.815.**

Every rung of the ladder is measured **from** B0 or B0+, so this is not a
constant offset — a contaminated floor makes `B0+ − B0`, `B1 − B0+` and
`B2 − B0+` uninterpretable.

**Fix** (`06_baselines.py`): add an explicit `log10(plasmid_length)` column read
from the reference FASTA, outside the pooled block; **B0 uses it alone**; keep
the current definition as a separate rung named **`B0_pooled`**; **assert** B0
has exactly one column and that its index is not part of the pooled vector — a
silent regression here is invisible in the output, which is how it survived the
first cohort.

**Deadline: before the v2 cohort is scored.**

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

- [ ] **A fixed and asserted** in `06_baselines.py`
- [ ] **Destructive cleanup done** (`rm -rf work/asm work/graphs`, etc.) — stale
      v1 assemblies would be labelled against v2 references, return all-absent,
      and enter the feature matrix as plausible data
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
- [ ] **Peak scratch stayed inside plan** — v1 was ~12 GB, v2 design stages
      ~50 GB before Phase B starts; confirm actual peak against `df` logs.

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

**14** merged into **8c**, same argument as 8. **6** subsumed by open item **B**.

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
