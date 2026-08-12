#!/bin/bash
# Is any read set truncated / malformed? A FASTQ must have a multiple of 4 lines
# and every record must have seq length == qual length.
set -uo pipefail
echo "wsl uptime(PID1 elapsed): $(ps -o etime= -p 1)"
echo
printf '%-18s %10s %12s %8s %s\n' TAG LINES BYTES "MOD4" STATUS
for S in "$HOME"/scratch/*/; do
  T=$(basename "$S")
  F="$S/reads.fastq"
  [ -s "$F" ] || { printf '%-18s %10s\n' "$T" "MISSING"; continue; }
  L=$(wc -l < "$F")
  B=$(stat -c %s "$F")
  M=$((L % 4))
  BAD=$(awk 'NR%4==2{s=length($0)} NR%4==0{if(length($0)!=s){print NR; exit}}' "$F")
  ST="ok"
  [ "$M" -ne 0 ] && ST="TRUNCATED(lines%4=$M)"
  [ -n "$BAD" ] && ST="$ST SEQ/QUAL_MISMATCH@line$BAD"
  printf '%-18s %10d %12d %8d %s\n' "$T" "$L" "$B" "$M" "$ST"
done
