# Plasmid recovery GNN

Predicting which plasmids a long-read assembler will fail to recover — and by
which mechanism — from the assembly graph alone, using GraphSAGE.

## Dataset setup

**The dataset is not in this repository** (see `.gitignore` — it is 7.1 GB). It
must be placed at `Dataset/` in the repository root before anything will run.
Every script hardcodes that path.

```
Dataset/
  sequences.fasta                              6.92 GB  raw NCBI plasmid query
  ecoli_10/ ecoli_10.zip  ecoli_ids.txt                 10 RefSeq E. coli genomes
  assemblies.tar.gz  assemblies/                        Wick reference assemblies (7 isolates)
  johnson-et-al-2023_supplementary-data-{1,2}.csv       Johnson recovery tables
  29584931/                                             NEKSUS accessions + contig summary
  datasets.exe                                          NCBI datasets CLI
```

Fastest route: **copy it from the shared drive.** Failing that, each component
is obtainable separately:

- **`ecoli_10/`** — exactly reproducible. The ten accessions are recorded below,
  so this survives even though `Dataset/` is untracked:
  `GCF_000597845.1 GCF_000599625.1 GCF_000599645.1 GCF_000599665.1
  GCF_000599685.1 GCF_000599705.1 GCF_000784925.1 GCF_000801165.1
  GCF_000801185.2 GCF_000814145.2`
  Fetch with the NCBI `datasets` CLI.
- **Johnson tables** — supplementary data 1 and 2 from Johnson et al. 2023.
- **NEKSUS** — figshare item `29584931`.
- **Wick** — the deposited `assemblies.tar.gz`.
- **`sequences.fasta`** — ⚠️ **the exact NCBI query string is not recorded
  anywhere.** It is an unfiltered nucleotide plasmid query yielding 72,556
  records / 7.3 Gbp, not PLSDB. Ask the maintainer for the query before trying
  to reconstruct it; a different query gives a different reference pool and the
  cohort will not reproduce. This is a known provenance gap — PROGRESS.md §12.6.

`01_filter_plasmids.py` does **not** download anything. It streams an existing
`sequences.fasta` and regenerates `work/refs/plasmid_manifest.csv`.

## Status and open items

- **`work/PROGRESS.md`** — current status: what is built, what was measured,
  what the numbers say, and the v2 cohort rebuild.
- **`work/PROPOSAL_FIXES.md`** — open items. Two are live, one of them blocking.
