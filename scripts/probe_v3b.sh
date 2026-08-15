#!/bin/bash
# v3 DESIGN PROBE ROUND 2 -- assemble ladder C under the two conditions most
# and least favourable to fragmentation.
#
#   r10_hq @ 30x  -- cleanest reads, most coverage. Directly comparable to
#                    round 1, so any change is attributable to ladder C's
#                    identity / copy-count / placement changes and not to the
#                    assembly conditions.
#   r9_raw  @ 15x -- ~91 % identity reads and half the depth: the condition in
#                    the cohort where fragmentation is most likely. If ladder C
#                    does not fragment HERE, it does not fragment anywhere.
#
# prep = rap (no --small_plasmid_bias) for both: every ladder-C carrier is
# 48-200 kb, far above the bias's range, so prep is irrelevant here and holding
# it fixed keeps the comparison clean.
#
# r9_raw needs ~9 GB and gets 6 threads -- it OOM-killed this VM three times
# during v1. Phase B is sequential. Do not parallelise.
set -uo pipefail
source /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/scripts/env.sh

DROOT=/mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work
PROBE=$DROOT/probe2
REF=$PROBE/sim/probe02_ref.fasta
S=$HOME/probe2_scratch
FRAGLEN=15000,13000
SEED=20260815

[ -s "$REF" ] || { echo "no reference -- run probe_v3b.py first"; exit 1; }
mkdir -p "$S" "$PROBE/asm"
echo "reference: $(grep -c '^>' "$REF") records  $(date +%H:%M:%S)"

# cell:  tag  arm  error_model  depth
CELLS="probe02_rap_d30:r10_hq:nanopore2023:30 probe02_rap_d15:r9_raw:nanopore2020:15"

# ---------- Phase A ----------------------------------------------------------
for C in $CELLS; do
  TAG=${C%%:*}; REST=${C#*:}; ARM=${REST%%:*}; REST=${REST#*:}
  EM=${REST%%:*}; DEPTH=${REST##*:}
  D="$S/$TAG"
  [ -s "$D/reads.fastq" ] && { echo "[$TAG] reads exist, skip"; continue; }
  mkdir -p "$D"
  ( badread simulate --reference "$REF" --quantity "${DEPTH}x" \
        --length "$FRAGLEN" --error_model "$EM" --qscore_model "$EM" \
        --seed "$SEED" > "$D/reads.fastq" 2> "$D/badread.err"
    echo "[$TAG] reads done $(du -h "$D/reads.fastq" | cut -f1) $(date +%H:%M:%S)" ) &
done
wait
echo "PHASE A COMPLETE $(date +%H:%M:%S)"

for C in $CELLS; do
  TAG=${C%%:*}
  awk 'NR%4==1{split($2,a,","); src=a[1]}
       NR%4==2{b[src]+=length($0); c[src]++}
       END{for(k in b) printf "%s\t%d\t%d\n", k, c[k], b[k]}' \
      "$S/$TAG/reads.fastq" | sort > "$S/$TAG/replicon_yield.tsv"
done

# ---------- Phase B: sequential ---------------------------------------------
for C in $CELLS; do
  TAG=${C%%:*}; REST=${C#*:}; ARM=${REST%%:*}
  D="$S/$TAG"; OUT="$PROBE/asm/$TAG"
  mkdir -p "$OUT"
  [ -s "$OUT/assembly_graph.gfa" ] && { echo "[$TAG] gfa exists, skip"; continue; }
  case "$ARM" in
    r9_raw) MODE="--nano-raw"; TH=6  ;;
    *)      MODE="--nano-hq";  TH=10 ;;
  esac
  rm -rf "$D/flye"
  t0=$(date +%s)
  flye $MODE "$D/reads.fastq" --out-dir "$D/flye" --threads $TH \
      > "$D/flye.stdout" 2> "$D/flye.stderr"
  rc=$?; t1=$(date +%s)
  if [ $rc -ne 0 ] || [ ! -s "$D/flye/assembly_graph.gfa" ]; then
    echo "[$TAG] FLYE FAILED rc=$rc after $((t1-t0))s"
    tail -8 "$D/flye.stderr" | sed 's/^/    /'
    continue
  fi
  cp "$D/flye/assembly_graph.gfa" "$D/flye/assembly_info.txt" \
     "$D/flye/assembly.fasta" "$D/replicon_yield.tsv" "$OUT/" 2>/dev/null
  SEG=$(grep -c '^S' "$OUT/assembly_graph.gfa")
  REAL=$(awk -F'\t' '$1=="L" && $2!=$4' "$OUT/assembly_graph.gfa" | wc -l)
  SELF=$(awk -F'\t' '$1=="L" && $2==$4' "$OUT/assembly_graph.gfa" | wc -l)
  echo "[$TAG] $MODE ok $((t1-t0))s  segments=$SEG  REAL-edges=$REAL  self-loops=$SELF"
done

echo "PROBE 2 COMPLETE $(date +%H:%M:%S)"
free -m | head -2
