#!/bin/bash
echo "--- os ---"
head -2 /etc/os-release
echo "--- python ---"
python3 --version
echo "--- tools ---"
for t in flye minimap2 badread seqkit rasusa conda mamba pip3 gcc make curl wget; do
  if command -v "$t" >/dev/null 2>&1; then echo "$t = $(command -v $t)"; else echo "$t = MISSING"; fi
done
echo "--- cpu ---"
nproc
echo "--- mem(GB) ---"
free -g | head -2
echo "--- disk ---"
df -h / 2>/dev/null | tail -1
echo "--- /mnt/d visible? ---"
ls /mnt/d/DNA-Sequencing/Dataset 2>&1 | head -8
