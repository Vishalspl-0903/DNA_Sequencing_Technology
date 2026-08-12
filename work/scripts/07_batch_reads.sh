#!/bin/bash
# Phase A -- generate reads for every (isolate, prep, depth) cell, in parallel.
#
# COHORT DESIGN
#   prep   lig = --small_plasmid_bias   (ligation-like: Badread does NOT
#                                        compensate for fragments longer than a
#                                        small circular plasmid, so small
#                                        plasmids are depleted size-dependently)
#          rap = no bias                (rapid-like: small plasmids retained)
#          Measured on sim01 at 5x: the 2.7 kb plasmid realises 27x relative
#          depth under rap and 1.9x under lig -- a 17x prep effect -- while the
#          85 kb plasmid is unaffected (3.0x vs 2.6x). This is the Wick
#          ligation/rapid contrast reproduced at zero bandwidth.
#
#   depth  30x and 15x (exact, set directly -- the simulation ladder is not a
#          subsample)
#
#   FRAGMENT LENGTH IS HELD FIXED across both preps, IDENTICALLY. Real ligation
#   and rapid preps also differ in read length, but confounding prep with N50
#   would make the prep effect uninterpretable. N50 is swept separately.
#
#   COHORT v2: FRAGLEN 6000,4000 -> 15000,13000.
#   At 6000,4000 the prep depletion peaked on sub-4 kb plasmids -- but those are
#   shorter than nearly every read, so Flye read the wrap-around as a repeat and
#   dropped them in ALL FOUR CELLS regardless of prep or depth. The 17x
#   depletion measured at the DEPTH level never reached the LABEL level: cell
#   failure rates came out 26.2/21.4/21.4/23.8%, flat, and only 5 of 84 pairs
#   flipped on prep. Depletion operates around the fragment length, so raising
#   it moves the affected band up onto 4-12 kb, where dropout is probabilistic
#   and where the reported ligation depletion actually operates.
#
#   CAVEAT, STATED NOT BURIED: raising the mean also raises wrap-around for the
#   4-12 kb band itself, which is the same mechanism that made sub-4 kb
#   deterministic at the old setting. The wide sd (13000) is what keeps the band
#   probabilistic -- fragments span roughly 2-28 kb, so a 4-12 kb plasmid gets a
#   mix of shorter and longer reads rather than uniformly longer ones. Whether
#   that holds is not assumed: the gate in 09_report.py requires cell failure
#   rates to SPREAD by >=10 pp, and 06_baselines.py refuses to analyse a cohort
#   with <20 discordant pairs. If the band goes deterministic again, both fire.
#
#   Badread error/qscore models are pinned per chemistry arm, never defaulted.
#
#   SEEDS: derived deterministically per (isolate, prep, depth) so every cell is
#   an independent draw AND the whole cohort reproduces from BASE_SEED. The seed
#   is written to the scratch dir and carried into the assembly manifest --
#   without that the read set cannot be regenerated after the reads are deleted.
#
# Badread is single-threaded; Flye is memory-bound and runs in Phase B.
set -uo pipefail

export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"
RUN="$MM run -n plasmid -r $MAMBA_ROOT_PREFIX"
DROOT="/mnt/d/DNA-Sequencing/work"
PAR="${PAR:-8}"
DEPTHS="${DEPTHS:-30 15}"
PREPS="${PREPS:-lig rap}"
FRAGLEN="${FRAGLEN:-15000,13000}"
BASE_SEED="${BASE_SEED:-20260811}"

mkdir -p "$HOME/scratch"

gen () {
  ISO="$1"; ARM="$2"; PREP="$3"; DEPTH="$4"; SEED="$5"
  TAG="${ISO}_${PREP}_d${DEPTH}"
  S="$HOME/scratch/$TAG"
  if [ -s "$S/reads.fastq" ]; then echo "[$TAG] reads exist, skip"; return 0; fi
  mkdir -p "$S"
  case "$ARM" in
    r10_hq) EM=nanopore2023 ;;
    r9_raw) EM=nanopore2020 ;;
    *)      EM=nanopore2023 ;;
  esac
  BIAS=""
  [ "$PREP" = "lig" ] && BIAS="--small_plasmid_bias"

  t0=$(date +%s)
  $RUN badread simulate \
      --reference "$DROOT/sim/${ISO}_ref.fasta" \
      --quantity "${DEPTH}x" \
      --length "$FRAGLEN" \
      --error_model "$EM" --qscore_model "$EM" \
      --seed "$SEED" $BIAS \
      > "$S/reads.fastq" 2> "$S/badread.raw"
  rc=$?
  t1=$(date +%s)
  tr '\r' '\n' < "$S/badread.raw" | grep -v '^Simulating:' | grep -v '^$' \
      | tail -30 > "$S/badread.log"
  rm -f "$S/badread.raw"

  # realised per-replicon coverage, straight from the read source tags
  # header looks like: @<uuid> chromosome,+strand,3269953-3306569 length=35806
  awk 'NR%4==1{split($2,a,","); src=a[1]}
       NR%4==2{bases[src]+=length($0); cnt[src]++}
       END{for(k in bases) printf "%s\t%d\t%d\n", k, cnt[k], bases[k]}' \
      "$S/reads.fastq" | sort > "$S/replicon_yield.tsv"

  printf '%s\n' "$ARM" > "$S/arm"
  printf '%s\n' "$PREP" > "$S/prep"
  printf '%s\n' "$ISO"  > "$S/isolate"
  printf '%s\n' "$FRAGLEN" > "$S/fraglen"
  # provenance: reads are deleted after graph extraction, so the seed and the
  # error model have to survive into the manifest or the cell is unreproducible
  printf '%s\n' "$SEED" > "$S/seed"
  printf '%s\n' "$EM"   > "$S/error_model"
  if [ $rc -ne 0 ]; then echo "[$TAG] BADREAD FAILED rc=$rc"; return 1; fi
  echo "[$TAG] ok $((t1-t0))s $(du -h "$S/reads.fastq" | cut -f1)"
}
export -f gen
export MM MAMBA_ROOT_PREFIX RUN DROOT FRAGLEN

i=0
# designs.json is {"provenance":..., "isolates":[...]} from cohort v2 and a bare
# list in v1; accept both so an old cohort can still be regenerated.
for line in $(python3 -c "
import json
d = json.load(open('$DROOT/sim/designs.json'))
for x in (d['isolates'] if isinstance(d, dict) else d):
    print(x['isolate']+','+x['arm'])
"); do
  ISO="${line%%,*}"; ARM="${line##*,}"
  i=$((i+1))
  for P in $PREPS; do
    for D in $DEPTHS; do
      # one distinct, reproducible seed per (isolate, prep, depth) cell.
      # i*100 leaves room for the prep offset (0/10/20) plus the depth, so no
      # two cells can collide; the whole cohort replays from BASE_SEED.
      case "$P" in lig) PO=0 ;; rap) PO=10 ;; *) PO=20 ;; esac
      CELL_SEED=$(( BASE_SEED + i * 100 + PO + D ))
      while [ "$(jobs -rp | wc -l)" -ge "$PAR" ]; do sleep 5; done
      gen "$ISO" "$ARM" "$P" "$D" "$CELL_SEED" &
    done
  done
done
wait
echo "PHASE A COMPLETE"
echo "read sets: $(ls -d "$HOME"/scratch/*/ 2>/dev/null | wc -l)"
du -sh "$HOME/scratch" | cut -f1
