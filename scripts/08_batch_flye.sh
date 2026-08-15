#!/bin/bash
# Phase B -- assemble every prepared read set with Flye.
#
# Flye mode is pinned by chemistry and never varied within a cohort:
#   r10_hq -> --nano-hq   (R10.4.1 / Dorado super-accuracy)
#   r9_raw -> --nano-raw  (R9.4.1  / Guppy)
#
# The env is activated by PATH (scripts/env.sh), NOT via `micromamba run`: the
# wrapper forwards signals to its process group, and `$(... | head -1)` in the
# manifest was SIGTERMing this script right after the first assembly finished.
#
# Reads are deleted after graph extraction (proposal Sec. II); the manifest
# retains provenance -- including the Badread seed and error model, without
# which a deleted read set cannot be regenerated.
#
# CONCURRENCY -- DO NOT PARALLELISE THIS LOOP.
# This loop is deliberately SEQUENTIAL: one Flye at a time, no `&`, no job pool.
# Flye is memory-bound, not core-bound. On R9.4.1-era reads (~91% identity,
# Badread nanopore2020) its k-mer index needs >=9 GB and it OOM-killed the whole
# WSL VM three times before the allocation was raised; even at 9 GB it peaked
# with 4 MB to spare. Peak also scales with --threads, which is why threads are
# pinned PER ARM (6 for r9_raw, 10 for r10_hq) rather than set globally.
# Cohort v2 is 200 assemblies, half of them on the R9 arm at ~190 s each.
# Concurrency on that arm is RAM / 9 GB -- on this machine that is ONE. Raising
# THREADS_R9 or backgrounding this loop reintroduces the OOM.
set -uo pipefail
source /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/scripts/env.sh

DROOT="/mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work"
KEEP_READS="${KEEP_READS:-0}"

FLYE_VERSION="$(flye --version 2>&1 | head -1)"

for S in "$HOME"/scratch/*/; do
  TAG=$(basename "$S")
  [ -s "$S/reads.fastq" ] || { echo "[$TAG] no reads, skip"; continue; }
  # completion marker written by gen() only after Badread exits
  if [ ! -s "$S/fraglen" ] || [ ! -s "$S/replicon_yield.tsv" ]; then
    echo "[$TAG] reads incomplete, skip"; continue
  fi
  OUT="$DROOT/asm/$TAG"
  if [ -s "$OUT/assembly_graph.gfa" ]; then echo "[$TAG] gfa exists, skip"; continue; fi

  ARM=$(cat "$S/arm" 2>/dev/null || echo r10_hq)
  PREP=$(cat "$S/prep" 2>/dev/null || echo lig)
  ISO=$(cat "$S/isolate" 2>/dev/null || echo "${TAG%%_*}")
  FRAG=$(cat "$S/fraglen" 2>/dev/null || echo "")
  CELLSEED=$(cat "$S/seed" 2>/dev/null || echo "null")
  EMODEL=$(cat "$S/error_model" 2>/dev/null || echo "")
  # Flye's peak memory scales with thread count, and the r9_raw arm
  # (nanopore2020, ~91% identity) is far hungrier than r10_hq: at 8.7 GB it
  # peaked with only 4 MB to spare at 8 threads / 15x. Give that arm fewer
  # threads for headroom; r10_hq is cheap and can have more.
  case "$ARM" in
    r9_raw) MODE="--nano-raw"; TH="${THREADS_R9:-6}"  ;;
    *)      MODE="--nano-hq";  TH="${THREADS_R10:-10}" ;;
  esac

  # never resume a partial work dir -- Flye would continue from inconsistent state
  rm -rf "$S/flye"
  mkdir -p "$OUT"

  t0=$(date +%s)
  flye $MODE "$S/reads.fastq" --out-dir "$S/flye" --threads "$TH" \
      > "$S/flye.stdout" 2> "$S/flye.stderr"
  rc=$?
  t1=$(date +%s)

  if [ $rc -ne 0 ] || [ ! -s "$S/flye/assembly_graph.gfa" ]; then
    echo "[$TAG] FLYE FAILED rc=$rc after $((t1-t0))s"
    tail -8 "$S/flye.stderr" | sed 's/^/    /'
    cp "$S/flye.stderr" "$OUT/flye.stderr" 2>/dev/null
    continue
  fi

  cp "$S/flye/assembly_graph.gfa" "$S/flye/assembly_info.txt" \
     "$S/flye/assembly.fasta" "$OUT/" 2>/dev/null
  cp "$S/badread.log" "$S/replicon_yield.tsv" "$OUT/" 2>/dev/null

  READBASES=$(awk 'NR%4==2{n+=length($0)} END{print n+0}' "$S/reads.fastq")
  DEPTH="${TAG##*_d}"
  BIASFLAG=false
  [ "$PREP" = "lig" ] && BIASFLAG=true
  NSEG=$(grep -c '^S' "$OUT/assembly_graph.gfa")

  cat > "$OUT/manifest.json" <<EOF
{"isolate":"$ISO","tag":"$TAG","arm":"$ARM","prep":"$PREP",
 "depth_fraction_x":$DEPTH,"fragment_length":"$FRAG",
 "badread_seed":$CELLSEED,"badread_error_model":"$EMODEL",
 "badread_qscore_model":"$EMODEL",
 "flye_mode":"$MODE","threads":$TH,"read_bases":${READBASES:-0},
 "seconds_flye":$((t1-t0)),"small_plasmid_bias":$BIASFLAG,
 "n_segments":$NSEG,"flye_version":"$FLYE_VERSION"}
EOF

  echo "[$TAG] ok $((t1-t0))s segments=$NSEG"
  [ "$KEEP_READS" = "0" ] && rm -rf "$S"
done

echo "PHASE B COMPLETE graphs=$(ls -d "$DROOT"/asm/*/assembly_graph.gfa 2>/dev/null | wc -l)"
