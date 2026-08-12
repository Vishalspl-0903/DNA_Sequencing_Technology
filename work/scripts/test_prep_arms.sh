#!/bin/bash
# Does the no-bias (rapid-prep-like) arm work for small plasmids, and does it
# actually change small-plasmid representation? Cheap low-depth probe.
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"
RUN="$MM run -n plasmid -r $MAMBA_ROOT_PREFIX"
REF=/mnt/d/DNA-Sequencing/work/sim/sim01_ref.fasta
T=$HOME/preptest; mkdir -p $T

probe () {
  NAME="$1"; shift
  echo "=== $NAME : badread $* ==="
  $RUN badread simulate --reference "$REF" --quantity 5x --seed 1 \
      --error_model nanopore2023 --qscore_model nanopore2023 "$@" \
      > "$T/$NAME.fastq" 2> "$T/$NAME.err"
  rc=$?
  if [ $rc -ne 0 ]; then
    echo "  FAILED rc=$rc"
    tr '\r' '\n' < "$T/$NAME.err" | grep -iv '^simulating' | tail -6 | sed 's/^/    /'
    return
  fi
  tr '\r' '\n' < "$T/$NAME.err" | grep -i 'depth' | head -6 | sed 's/^/    /'
  echo "  bases: $(awk 'NR%4==2{n+=length($0)} END{print n}' $T/$NAME.fastq)"
}

probe rapid_default                                   # no bias, default 15k frags
probe rapid_short   --length 6000,4000                # no bias, shorter frags
probe ligation_short --length 6000,4000 --small_plasmid_bias
