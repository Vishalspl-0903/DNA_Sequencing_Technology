#!/bin/bash
# Assemble ONE read set in the foreground while sampling memory, to find out
# whether Flye's peak is what takes the WSL VM down.
set -uo pipefail
source /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/scripts/env.sh

TH="${TH:-8}"
S=$(ls -d "$HOME"/scratch/*/ 2>/dev/null | head -1)
TAG=$(basename "$S")
echo "assembling $TAG with $TH threads"
echo "total mem: $(free -m | awk 'NR==2{print $2}') MB"

rm -rf "$S/flye"
( while true; do
    free -m | awk 'NR==2{printf "  mem used=%s free=%s avail=%s\n",$3,$4,$7}'
    sleep 15
  done ) &
SAMPLER=$!

t0=$(date +%s)
flye --nano-hq "$S/reads.fastq" --out-dir "$S/flye" --threads "$TH" \
    > "$S/flye.stdout" 2> "$S/flye.stderr"
rc=$?
t1=$(date +%s)
kill $SAMPLER 2>/dev/null

echo "rc=$rc in $((t1-t0))s"
if [ -s "$S/flye/assembly_graph.gfa" ]; then
  echo "OK segments=$(grep -c '^S' "$S/flye/assembly_graph.gfa")"
else
  echo "NO GFA"; tail -5 "$S/flye.stderr" | sed 's/^/    /'
fi
echo "peak-ish free after: $(free -m | awk 'NR==2{print $7}') MB available"
