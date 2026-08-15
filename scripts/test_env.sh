#!/bin/bash
source /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/scripts/env.sh
echo "--- tools resolve on PATH ---"
for t in flye minimap2 seqkit rasusa badread python; do
  printf '%-10s %s\n' "$t" "$(command -v $t || echo MISSING)"
done
echo "--- the exact construct that was killing the driver ---"
V="$(flye --version 2>&1 | head -1)"
echo "flye --version | head -1  => $V"
echo "SURVIVED the pipe-to-head"
echo "--- still alive after a second one ---"
W="$(minimap2 --version 2>&1 | head -1)"
echo "minimap2 => $W"
echo "TEST COMPLETE"
