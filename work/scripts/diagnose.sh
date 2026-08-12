#!/bin/bash
echo "--- running procs ---"
ps -eo pid,ppid,pcpu,rss,etime,comm --sort=-rss | head -12
echo "--- memory ---"
free -m | head -3
echo "--- OOM / kill evidence ---"
dmesg 2>/dev/null | tail -25 | grep -iE "oom|kill|memory" || echo "(dmesg unavailable or no oom lines)"
echo "--- scratch remaining ---"
ls -d "$HOME"/scratch/*/ 2>/dev/null | wc -l
echo "--- graphs done ---"
ls -d /mnt/d/DNA-Sequencing/work/asm/*/assembly_graph.gfa 2>/dev/null | wc -l
echo "--- partial flye dirs ---"
for d in "$HOME"/scratch/*/flye; do [ -d "$d" ] && echo "  $(basename "$(dirname "$d")")"; done
echo "--- disk ---"
df -h / | tail -1
