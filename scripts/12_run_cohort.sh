#!/bin/bash
# Cohort driver -- Phase A and Phase B interleaved in isolate slices.
#
# WHY THIS EXISTS
# 07_batch_reads.sh stages every read set it generates, and 08_batch_flye.sh
# only deletes a read set after its graph is extracted. Running A to completion
# before B therefore peaks at the FULL cohort's read volume. Measured from
# cohort v1: 10 GB of FASTQ for 12 isolates, so 50 isolates is ~42 GB. That
# scratch lives in the WSL ext4.vhdx on C:, which has 40.7 GB free, and the
# vhdx does not shrink again once grown. It does not fit.
#
# Slicing A->B->A->B in blocks of SLICE isolates holds the peak at
# SLICE * ~0.85 GB. At the default 10 that is ~8.5 GB.
#
# REPRODUCIBILITY IS NOT AFFECTED. 07 derives each cell's seed from the
# isolate's position in the full designs.json list and advances that counter
# for every isolate whether or not it falls in the slice, so a sliced run and
# an unsliced run generate byte-identical reads.
#
# CONCURRENCY. Phase A and Phase B never overlap. Badread is light but Flye on
# the R9 arm needs ~9 GB of a 9 GB VM and OOM-killed this VM three times during
# v1; nothing else runs while it does.
set -uo pipefail

DROOT="/mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work"

# Positional args as well as env vars: wsl.exe does not reliably preserve a
# quoted `bash -c '...'` string across the interop boundary, so `START=21 wsl
# ... 12_run_cohort.sh` silently loses the assignment. Resuming after an
# interrupted slice must be `bash 12_run_cohort.sh 21`.
#   $1 = first isolate (1-based, default 1)
#   $2 = isolates per slice (default 10)
#   $3 = Badread parallelism in Phase A (default 8)
#
# $3 only affects Phase A. Badread is single-threaded and light on memory, and
# nothing else runs during Phase A, so on 24 cores this is safe to raise --
# measured slice 1 took 58 min at PAR=8. Phase B is NOT affected and must not
# be: Flye needs ~9 GB of a 9 GB VM and is sequential by design.
START="${1:-${START:-1}}"
SLICE="${2:-${SLICE:-10}}"
export PAR="${3:-${PAR:-8}}"

N=$(python3 -c "
import json
d = json.load(open('$DROOT/sim/designs.json'))
print(len(d['isolates'] if isinstance(d, dict) else d))
")
[ -z "$N" ] && { echo "cannot read designs.json -- run 02_build_sim_refs.py first"; exit 1; }

echo "=============================================================="
echo "COHORT DRIVER  isolates=$N  slice=$SLICE  starting at $START"
echo "=============================================================="

t_start=$(date +%s)
from=$START
while [ "$from" -le "$N" ]; do
  to=$(( from + SLICE - 1 ))
  [ "$to" -gt "$N" ] && to=$N

  echo
  echo "########## SLICE $from-$to : PHASE A  $(date +%H:%M:%S) ##########"
  ISO_FROM=$from ISO_TO=$to bash "$DROOT/scripts/07_batch_reads.sh"

  echo
  echo "########## SLICE $from-$to : PHASE B  $(date +%H:%M:%S) ##########"
  bash "$DROOT/scripts/11_assemble_all.sh"

  graphs=$(ls -d "$DROOT"/asm/*/assembly_graph.gfa 2>/dev/null | wc -l)
  left=$(ls -d "$HOME"/scratch/*/ 2>/dev/null | wc -l)
  free_c=$(df -BG --output=avail /  | tail -1 | tr -d ' G')
  echo "--- slice $from-$to done: graphs=$graphs  read_sets_left=$left"
  echo "    vm free disk=${free_c}G  elapsed=$(( ($(date +%s) - t_start) / 60 ))m"

  # a read set that survives Phase B means its assembly failed; stop rather
  # than pile the next slice's 8.5 GB on top of it
  if [ "$left" -ne 0 ]; then
    echo "!! $left read set(s) unassembled after slice $from-$to -- STOPPING."
    echo "   inspect, then resume with: START=$from bash $0"
    exit 1
  fi

  from=$(( to + 1 ))
done

echo
echo "=============================================================="
echo "COHORT COMPLETE  graphs=$(ls -d "$DROOT"/asm/*/assembly_graph.gfa 2>/dev/null | wc -l)/$(( N * 4 ))"
echo "elapsed $(( ($(date +%s) - t_start) / 60 )) min   $(date +%H:%M:%S)"
echo "=============================================================="
