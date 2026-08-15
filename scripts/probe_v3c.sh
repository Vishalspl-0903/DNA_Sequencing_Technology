#!/bin/bash
# v3 DESIGN PROBE ROUND 3 -- assemble the length x copy-number map over the
# FULL 2x2 the cell-spread gate is computed on, plus one r9_raw check.
#
#   lig/rap x 15x/30x on r10_hq  -- the four cells the gate averages over.
#     The gate failed at 5.5 pp in v2 and this is the only way to predict it
#     before committing another 12 h.
#
#   rap @ 15x on r9_raw          -- the cohort is half r9_raw. The ligation
#     bias is a Badread read-SAMPLING effect and is chemistry-independent, but
#     the error model changes how hard assembly is, which can move the window.
#     One cell is enough to see whether it moves a lot or a little.
#
# Phase B is sequential: r9_raw needs ~9 GB of a 9 GB VM.
set -uo pipefail
source /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/scripts/env.sh

DROOT=/mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work
PROBE=$DROOT/probe3
REF=$PROBE/sim/probe03_ref.fasta
S=$HOME/probe3_scratch
FRAGLEN=15000,13000
SEED=20260816

[ -s "$REF" ] || { echo "no reference -- run probe_v3c.py first"; exit 1; }
mkdir -p "$S" "$PROBE/asm"
echo "reference: $(grep -c '^>' "$REF") records  $(date +%H:%M:%S)"

# tag : prep : depth : arm : error_model
CELLS="probe03_lig_d30:lig:30:r10_hq:nanopore2023
probe03_rap_d30:rap:30:r10_hq:nanopore2023
probe03_lig_d15:lig:15:r10_hq:nanopore2023
probe03_rap_d15:rap:15:r10_hq:nanopore2023
probe03r9_rap_d15:rap:15:r9_raw:nanopore2020"

# ---------- Phase A: all cells in parallel ----------------------------------
for C in $CELLS; do
  IFS=: read -r TAG PREP DEPTH ARM EM <<< "$C"
  D="$S/$TAG"
  [ -s "$D/reads.fastq" ] && { echo "[$TAG] reads exist, skip"; continue; }
  mkdir -p "$D"
  BIAS=""
  [ "$PREP" = "lig" ] && BIAS="--small_plasmid_bias"
  ( badread simulate --reference "$REF" --quantity "${DEPTH}x" \
        --length "$FRAGLEN" --error_model "$EM" --qscore_model "$EM" \
        --seed "$SEED" $BIAS > "$D/reads.fastq" 2> "$D/badread.err"
    echo "[$TAG] reads done $(du -h "$D/reads.fastq" | cut -f1) $(date +%H:%M:%S)" ) &
done
wait
echo "PHASE A COMPLETE $(date +%H:%M:%S)"

for C in $CELLS; do
  IFS=: read -r TAG _ _ _ _ <<< "$C"
  awk 'NR%4==1{split($2,a,","); src=a[1]}
       NR%4==2{b[src]+=length($0); c[src]++}
       END{for(k in b) printf "%s\t%d\t%d\n", k, c[k], b[k]}' \
      "$S/$TAG/reads.fastq" | sort > "$S/$TAG/replicon_yield.tsv"
done

# ---------- Phase B: sequential ---------------------------------------------
for C in $CELLS; do
  IFS=: read -r TAG PREP DEPTH ARM EM <<< "$C"
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
  echo "[$TAG] $MODE ok $((t1-t0))s segments=$(grep -c '^S' "$OUT/assembly_graph.gfa")"
done

echo "PROBE 3 COMPLETE $(date +%H:%M:%S)"
du -sh "$S"
