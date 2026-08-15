#!/bin/bash
# Cohort v2 smoke test -- isolates 2 and 3 only, end to end.
#
# Deliberately covers both chemistry arms in one pass:
#   sim02  r9_raw  --nano-raw  (the arm that OOM-killed the VM three times)
#   sim03  r10_hq  --nano-hq   (carries frag-role IS repeats)
#
# Invoked as a FILE, not `bash -c '...'`: wsl.exe does not reliably preserve
# a quoted -c string across the interop boundary (the shell variable inside it
# came back empty twice), so every WSL entry point here is an absolute path to
# a real script.
set -uo pipefail
R=/mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work

export ISO_FROM=2
export ISO_TO=3

echo "########## PHASE A  $(date +%H:%M:%S) ##########"
bash "$R/scripts/07_batch_reads.sh" || { echo "PHASE A FAILED"; exit 1; }

echo
echo "########## seeds written (reproducibility check) ##########"
for d in "$HOME"/scratch/*/; do
  printf "  %-20s seed=%s arm=%s model=%s frag=%s\n" \
    "$(basename "$d")" "$(cat "$d/seed")" "$(cat "$d/arm")" \
    "$(cat "$d/error_model")" "$(cat "$d/fraglen")"
done

echo
echo "########## PHASE B  $(date +%H:%M:%S) ##########"
bash "$R/scripts/11_assemble_all.sh"

echo
echo "########## RESULT ##########"
echo "graphs: $(ls -d "$R"/asm/*/assembly_graph.gfa 2>/dev/null | wc -l)"
echo "read sets left: $(ls -d "$HOME"/scratch/*/ 2>/dev/null | wc -l)"
free -m | head -2
df -h / | tail -1
