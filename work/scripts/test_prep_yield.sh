#!/bin/bash
# Measure realised per-replicon coverage in each prep arm.
T=$HOME/preptest
echo "head of a read name:"; head -1 $T/rapid_default.fastq
for f in rapid_default rapid_short ligation_short; do
  echo "=== $f ==="
  awk 'NR%4==1{split($2,a,"+"); split(a[1],b,"-"); src=b[1]; sub(/,.*/,"",src)}
       NR%4==2{bases[src]+=length($0); n[src]++}
       END{for(s in bases) printf "  %-14s reads=%-7d bases=%-12d\n", s, n[s], bases[s]}' \
      "$T/$f.fastq" | sort
done
