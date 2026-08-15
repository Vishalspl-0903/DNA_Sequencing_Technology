#!/bin/bash
# v3 DESIGN PROBE -- assemble the probe isolate under both preps.
#
# Two cells only: lig and rap, both at 30x, both r10_hq. Deliberate choices:
#
#   30x not 15x -- low depth causes fragmentation through coverage gaps, which
#   would confound the repeat-length question this probe exists to answer.
#   30x is the clean test: anything that fragments here fragments because of
#   the repeat, not because Flye ran out of reads.
#
#   r10_hq not r9_raw -- cleaner reads RESOLVE repeats more easily, so r10_hq
#   is the conservative case for ladder B. A repeat that defeats --nano-hq also
#   defeats --nano-raw. It is also 2-3x faster and needs half the memory.
#
#   FRAGLEN 15000,13000 -- identical to the cohort. The whole point is to
#   measure the ladders AT THE COHORT'S READ LENGTH.
#
# Writes to work/probe/, never work/asm -- the v2 corpus is untouched.
set -uo pipefail
source /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/scripts/env.sh

DROOT=/mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work
PROBE=$DROOT/probe
REF=$PROBE/sim/probe01_ref.fasta
S=$HOME/probe_scratch
FRAGLEN=15000,13000
DEPTH=30
SEED=20260814

[ -s "$REF" ] || { echo "no reference -- run probe_v3.py first"; exit 1; }
mkdir -p "$S" "$PROBE/asm"

echo "reference: $(grep -c '^>' "$REF") records, $(date +%H:%M:%S)"

# ---------- Phase A: both preps in parallel (badread is single-threaded) -----
for PREP in lig rap; do
  D="$S/probe01_${PREP}_d${DEPTH}"
  if [ -s "$D/reads.fastq" ]; then echo "[$PREP] reads exist, skip"; continue; fi
  mkdir -p "$D"
  BIAS=""
  [ "$PREP" = "lig" ] && BIAS="--small_plasmid_bias"
  ( badread simulate --reference "$REF" --quantity "${DEPTH}x" \
        --length "$FRAGLEN" --error_model nanopore2023 --qscore_model nanopore2023 \
        --seed "$SEED" $BIAS > "$D/reads.fastq" 2> "$D/badread.err"
    echo "[$PREP] reads done $(du -h "$D/reads.fastq" | cut -f1) $(date +%H:%M:%S)" ) &
done
wait
echo "PHASE A COMPLETE $(date +%H:%M:%S)"

# realised per-replicon yield, straight from the read source tags -- this is
# the DEPTH-level prep effect, before assembly gets a say
for PREP in lig rap; do
  D="$S/probe01_${PREP}_d${DEPTH}"
  awk 'NR%4==1{split($2,a,","); src=a[1]}
       NR%4==2{b[src]+=length($0); c[src]++}
       END{for(k in b) printf "%s\t%d\t%d\n", k, c[k], b[k]}' \
      "$D/reads.fastq" | sort > "$D/replicon_yield.tsv"
done

# ---------- Phase B: sequential, one Flye at a time --------------------------
for PREP in lig rap; do
  D="$S/probe01_${PREP}_d${DEPTH}"
  OUT="$PROBE/asm/probe01_${PREP}_d${DEPTH}"
  mkdir -p "$OUT"
  if [ -s "$OUT/assembly_graph.gfa" ]; then echo "[$PREP] gfa exists, skip"; continue; fi
  rm -rf "$D/flye"
  t0=$(date +%s)
  flye --nano-hq "$D/reads.fastq" --out-dir "$D/flye" --threads 10 \
      > "$D/flye.stdout" 2> "$D/flye.stderr"
  rc=$?
  t1=$(date +%s)
  if [ $rc -ne 0 ] || [ ! -s "$D/flye/assembly_graph.gfa" ]; then
    echo "[$PREP] FLYE FAILED rc=$rc after $((t1-t0))s"
    tail -8 "$D/flye.stderr" | sed 's/^/    /'
    continue
  fi
  cp "$D/flye/assembly_graph.gfa" "$D/flye/assembly_info.txt" \
     "$D/flye/assembly.fasta" "$D/replicon_yield.tsv" "$OUT/" 2>/dev/null
  echo "[$PREP] flye ok $((t1-t0))s segments=$(grep -c '^S' "$OUT/assembly_graph.gfa") links=$(grep -c '^L' "$OUT/assembly_graph.gfa")"
done

echo "PROBE COMPLETE $(date +%H:%M:%S)"
du -sh "$S"
