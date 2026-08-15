#!/bin/bash
echo "--- read sets in scratch ---"
for d in "$HOME"/scratch/*/; do
  t=$(basename "$d")
  if [ -s "$d/reads.fastq" ]; then
    printf '%-16s %8s\n' "$t" "$(du -h "$d/reads.fastq" | cut -f1)"
  else
    printf '%-16s %8s\n' "$t" "(pending)"
  fi
done
echo "--- assembled graphs on /mnt/d ---"
ls -d /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/asm/*/ 2>/dev/null | wc -l
for g in /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/asm/*/assembly_graph.gfa; do
  [ -s "$g" ] || continue
  printf '%-16s segments=%s\n' "$(basename "$(dirname "$g")")" "$(grep -c '^S' "$g")"
done 2>/dev/null
echo "--- running procs ---"
ps -eo pcpu,rss,etime,comm --sort=-pcpu | head -12
echo "--- memory (MB) ---"
free -m | head -2
