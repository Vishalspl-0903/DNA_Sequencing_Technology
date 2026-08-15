#!/bin/bash
# Run Flye writing its log straight to /mnt/d so the evidence survives the VM
# going down. Also drop a heartbeat file so we can see how far it got.
set -uo pipefail
source /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/scripts/env.sh

TAG="${TAG:-sim02_lig_d15}"
TH="${TH:-4}"
S="$HOME/scratch/$TAG"
# use the mode actually pinned for this isolate's chemistry, not a hardcoded one
ARM=$(cat "$S/arm" 2>/dev/null || echo r10_hq)
case "$ARM" in
  r9_raw) MODE="--nano-raw" ;;
  *)      MODE="--nano-hq"  ;;
esac
echo "arm=$ARM mode=$MODE"
D=/mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/logs/flyetest
mkdir -p "$D"
rm -rf "$S/flye"

echo "start $(date -u +%H:%M:%S) tag=$TAG threads=$TH mode=$MODE" > "$D/heartbeat.txt"
# sample fast (1s) -- a 10s sampler missed the spike that killed the VM
( while true; do
    echo "$(date -u +%H:%M:%S) mem_avail=$(free -m | awk 'NR==2{print $7}')" >> "$D/heartbeat.txt"
    sleep 1
  done ) &
HB=$!

flye $MODE "$S/reads.fastq" --out-dir "$S/flye" --threads "$TH" \
     > "$D/flye.stdout" 2> "$D/flye.stderr"
rc=$?
kill $HB 2>/dev/null
echo "FLYE_RC=$rc $(date -u +%H:%M:%S)" >> "$D/heartbeat.txt"
if [ -s "$S/flye/assembly_graph.gfa" ]; then
  echo "GFA_OK segments=$(grep -c '^S' "$S/flye/assembly_graph.gfa")" >> "$D/heartbeat.txt"
fi
echo "done rc=$rc"
