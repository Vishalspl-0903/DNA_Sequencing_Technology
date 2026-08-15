# Cohort v1 — archived corpus

**Archived:** 2026-08-13 · **Superseded by:** cohort v2 (`02_build_sim_refs.py`,
`COHORT_VERSION = 2`)

This is the complete 12-isolate cohort described in `PROGRESS.md` §4, §7, §8b
and §9. It was moved here — not deleted — because the assemblies are ~7 hours of
Flye and cannot be regenerated without regenerating the reads, which were
deleted after graph extraction (proposal Sec. II).

## Why it is here and not in `work/asm`

Nothing downstream is keyed to a cohort version. `04_label.py` and
`05_features.py` walk *every* directory under `work/asm/` and label it against
`work/sim/<isolate>_ref.fasta`. After the v2 rebuild, `sim01_ref.fasta` is a
**different isolate with different plasmids**, so a surviving v1
`asm/sim01_lig_d15` would be labelled against plasmids it never contained, come
back all-`absent`, and flow into `node_features.csv` and the baseline ladder as
if it were real data.

Both discovery calls are `os.listdir(ROOT/"asm")` and `os.listdir(ROOT/"graphs")`
— one level, no recursion — so `cohort_v1/` is invisible to the pipeline.
`check_stale_cohort()` in `02_build_sim_refs.py` likewise inspects only
`work/asm`, `work/graphs` and `work/sim/*_ref.fasta`, so this directory does not
trip the rebuild guard.

## Contents

| path | what | n |
|---|---|---|
| `asm/<tag>/` | `assembly_graph.gfa`, `assembly_info.txt`, `manifest.json`, `replicon_yield.tsv` | 49 tags |
| `graphs/<tag>/` | `plasmid_labels.csv`, `node_labels.csv`, `aln.paf` | 49 tags |
| `sim/` | `sim01..sim12_ref.fasta`, `ground_truth.csv`, `designs.json` | 12 isolates |
| `results/` | `node_features.csv`, `graph_diagnostics.json`, `pyg_dataset.pt`, `label_validation.csv`, `baselines.json` | — |

49 tags = 12 isolates × {lig, rap} × {30×, 15×} = 48 cohort cells, plus the
`sim01_d50` pilot. Total 786 MB.

## Provenance of these exact bytes

Copied from `D:\DNA-Sequencing\work`, which is the original working directory
Flye wrote into. **That copy, not the one previously in `work/asm`, is the
authoritative one.** The repository was checked out with `core.autocrlf=true`
and without a `.gitattributes`, so the working-tree copies of `*.gfa`, `*.fasta`
and `*.paf` had been rewritten with CRLF line endings. Every archived GFA was
verified to match its CRLF-stripped counterpart exactly (49/49, zero
mismatches), and the archive contains no CR bytes. `.gitattributes` now marks
these extensions `binary` so a future checkout cannot repeat it.

`results/baselines.json` here is the **real v1 ladder** (§9.4). The
`work/results/baselines.json` left in place is the v2 adequacy-gate stub written
when `06_baselines.py` refused to analyse this cohort — see §11.

## Status — do not quote these numbers

Cohort v1 **fails all four adequacy gates** (`PROGRESS.md` §11):

| gate | floor | v1 |
|---|---|---|
| discordant plasmids (prep or depth) | ≥ 20 | 7 |
| fragmented units, interventional | ≥ 20 | 8 |
| absorbed units, interventional | ≥ 20 | 1 |
| cell failure-rate spread | ≥ 10 pp | 6.5 pp |

The baseline ladder computed on it (B0 0.853 / B0+ 0.930 / B1 0.975 / B2 0.942)
is pipeline validation, not evidence. It is retained as the record of *why* the
rebuild happened — an underpowered cohort found inadequate by criteria fixed in
advance — and for regression-checking the labeller and feature extractor against
known outputs.
