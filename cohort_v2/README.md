# Cohort v2 — archived corpus

**Built:** 2026-08-13/14 (11.6 h of Flye) · **Archived:** 2026-08-14 ·
**Superseded by:** cohort v3

50 isolates × {lig, rap} × {30×, 15×} = **200 assemblies, zero Flye failures**,
696 plasmid label units. This is the corpus described in `PROGRESS.md` §10–§11
and it is kept because it is the evidence for every design decision in v3.

## Why it was superseded

It failed 2 of its 4 pre-registered gates, and the diagnosis is the useful part.

| gate | floor | v2 | |
|---|---|---|---|
| discordant plasmids | ≥ 20 | **27** | pass |
| absorbed units | ≥ 20 | **40** | pass |
| fragmented units | ≥ 20 | 19 | fail |
| cell failure spread | ≥ 10 pp | 5.5 pp | fail |

**Cell spread.** The `small` stratum was 4–12 kb × copy 8–30, which is three
regimes glued together: below ~8.4 kb units fail in both preps, above ~10.9 kb
they pass in both, and only the middle flips. Averaging them gave a number prep
could not move. v3 sets the band from measurement (`work/probe3`).

**Fragmented.** The induced mechanism converted **0 / 148**. Three probes then
established it is not inducible at all in this simulator — repeat length 2–25 kb,
2–6 copies in plasmid and 2–4 in chromosome, identity 0.990 and 1.000, r10_hq
and r9_raw, 30× and 15×: **0 / 164**. Nor is depth the lever — medium/large
units below 1.5× relative depth came out 56 recovered vs 6 fragmented, and
fragmented sits at a *higher* median relative depth than recovered (3.99 vs
2.87). The `fragmented` floor was therefore **withdrawn**, not relaxed.

## The finding that outlives the cohort

**45 real edges across 200 graphs; 175 of 200 edgeless; 698 of 768 nodes
isolated (91 %).** Three probes built specifically to create branch points
produced one edge between them. Badread samples uniformly from a perfect
circular template and Flye traverses that circle correctly whatever is inside
it, so the simulation arm does not generate assembly-graph topology.

A message-passing model on this corpus is an MLP with global pooling. The
relational-topology hypothesis (B3/B4) cannot be tested on simulated data — it
needs real assembly graphs. That is a result, and it is why effort moved to the
real cohorts.

## What v2 does support

- the end-to-end pipeline, on 200 assemblies with zero failures
- a large, real prep effect at the depth level (~10× ligation depletion)
- `absorbed` as an inducible mechanism: 36/40 = **90 % conversion**
- 27 discordant matched pairs, the design's primary contribution
- Arm C (fragmented/absorbed only) with **59 positives**, against v1's 9

## Contents

| path | | n |
|---|---|---|
| `asm/<tag>/` | GFA, `assembly_info.txt`, `manifest.json`, `replicon_yield.tsv` | 200 |
| `graphs/<tag>/` | `plasmid_labels.csv`, `node_labels.csv`, `aln.paf` | 200 |
| `sim/` | `sim01..sim50_ref.fasta`, `ground_truth.csv`, `designs.json` | 50 |
| `results/` | node features, diagnostics, PyG dataset, label validation, baselines | — |

3.2 GB. Sits outside `asm/` and `graphs/`, both discovered by a non-recursive
`listdir`, so it is invisible to the pipeline and does not trip
`check_stale_cohort()`. See also `cohort_v1/` for the 12-isolate predecessor.

**`baselines.json` here is the v2 gate-failure stub, not a ladder** — the cohort
never passed its gates, so `06_baselines.py` correctly refused to score it.
