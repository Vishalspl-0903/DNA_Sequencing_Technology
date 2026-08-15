#!/bin/bash
# Flye's disjointig-assembly stage is the memory hog, and on the r9_raw arm
# (nanopore2020, ~91% identity) it exhausts all 6.8 GB and OOMs the WSL VM.
# --asm-coverage uses only the longest N x of reads for that stage, which is the
# documented memory reduction; it requires --genome-size.
set -uo pipefail
source /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/scripts/env.sh

TAG="${TAG:-sim02_lig_d15}"
TH="${TH:-8}"
ASMCOV="${ASMCOV:-20}"
GSIZE="${GSIZE:-5m}"
S="$HOME/scratch/$TAG"
ARM=$(cat "$S/arm" 2>/dev/null || echo r10_hq)
case "$ARM" in
  r9_raw) MODE="--nano-raw" ;;
  *)      MODE="--nano-hq"  ;;
esac
D=/mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/logs/flyetest
mkdir -p "$D"
rm -rf "$S/flye"

echo "tag=$TAG arm=$ARM mode=$MODE threads=$TH asm-coverage=$ASMCOV"
echo "start $(date -u +%H:%M:%S) $TAG $MODE asmcov=$ASMCOV" > "$D/lowmem_hb.txt"
( while true; do
    echo "$(date -u +%H:%M:%S) avail=$(free -m | awk 'NR==2{print $7}')" >> "$D/lowmem_hb.txt"
    sleep 2
  done ) &
HB=$!

COVARGS=()
if [ "$ASMCOV" != "0" ]; then
  COVARGS=(--genome-size "$GSIZE" --asm-coverage "$ASMCOV")
fi

t0=$(date +%s)
flye $MODE "$S/reads.fastq" --out-dir "$S/flye" --threads "$TH" \
     "${COVARGS[@]}" \
     > "$D/lowmem.stdout" 2> "$D/lowmem.stderr"
rc=$?
t1=$(date +%s)
kill $HB 2>/dev/null
echo "RC=$rc secs=$((t1-t0))" >> "$D/lowmem_hb.txt"
echo "rc=$rc in $((t1-t0))s"
if [ -s "$S/flye/assembly_graph.gfa" ]; then
  echo "GFA OK segments=$(grep -c '^S' "$S/flye/assembly_graph.gfa")"
else
  echo "NO GFA"; tail -6 "$D/lowmem.stderr" | sed 's/^/    /'
fi
