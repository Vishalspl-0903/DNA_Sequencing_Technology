#!/bin/bash
# Phase 0 / step 3 -- one simulated isolate: reference -> reads -> Flye -> GFA.
#
# Usage: 03_sim_assemble.sh <isolate> <arm> <depth_x> [seed]
#   arm = r10_hq  -> Badread nanopore2023 + Flye --nano-hq   (R10.4.1 / Dorado sup)
#   arm = r9_raw  -> Badread nanopore2020 + Flye --nano-raw  (R9.4.1  / Guppy)
#
# Badread models are pinned per arm, never defaulted (proposal Sec. II): the
# default is nanopore2023, which would silently make the R9.4.1 arm R10.4.1-like
# and erase the chemistry contrast.
#
# --small_plasmid_bias is ON deliberately. Without it Badread inflates the depth
# of small circular plasmids to compensate for fragment-length misses, which
# erases the very dropout we are trying to predict. With it, small plasmids are
# depleted exactly as they are in a real ligation prep.
#
# Work happens on the WSL native filesystem (fast); only the graph and the
# assembly metadata are copied back to /mnt/d. Reads are deleted after graph
# extraction, per the proposal's manifest policy.
set -euo pipefail

ISO="$1"; ARM="$2"; DEPTH="$3"; SEED="${4:-1}"
THREADS="${THREADS:-16}"

export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"
RUN="$MM run -n plasmid -r $MAMBA_ROOT_PREFIX"

DROOT="/mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work"
REF="$DROOT/sim/${ISO}_ref.fasta"
TAG="${ISO}_d${DEPTH}"
SCRATCH="$HOME/scratch/$TAG"
OUTDIR="$DROOT/asm/$TAG"

case "$ARM" in
  r10_hq) EMODEL=nanopore2023; QMODEL=nanopore2023; FLYEMODE="--nano-hq"  ;;
  r9_raw) EMODEL=nanopore2020; QMODEL=nanopore2020; FLYEMODE="--nano-raw" ;;
  *) echo "unknown arm: $ARM" >&2; exit 2 ;;
esac

mkdir -p "$SCRATCH" "$OUTDIR"
echo "[$TAG] arm=$ARM depth=${DEPTH}x model=$EMODEL flye=$FLYEMODE threads=$THREADS"

t0=$(date +%s)
echo "[$TAG] simulating reads..."
$RUN badread simulate \
    --reference "$REF" \
    --quantity "${DEPTH}x" \
    --error_model "$EMODEL" \
    --qscore_model "$QMODEL" \
    --seed "$SEED" \
    --small_plasmid_bias \
    > "$SCRATCH/reads.fastq" 2> "$SCRATCH/badread.log"
t1=$(date +%s)
READBASES=$(awk 'NR%4==2{n+=length($0)} END{print n}' "$SCRATCH/reads.fastq")
echo "[$TAG] reads done in $((t1-t0))s, ${READBASES} bp"

echo "[$TAG] assembling with Flye..."
$RUN flye $FLYEMODE "$SCRATCH/reads.fastq" \
    --out-dir "$SCRATCH/flye" \
    --threads "$THREADS" \
    > "$SCRATCH/flye.stdout" 2> "$SCRATCH/flye.stderr" || {
      echo "[$TAG] FLYE FAILED"; tail -20 "$SCRATCH/flye.stderr"; exit 1; }
t2=$(date +%s)
echo "[$TAG] flye done in $((t2-t1))s"

cp "$SCRATCH/flye/assembly_graph.gfa"  "$OUTDIR/" 2>/dev/null || echo "no gfa!"
cp "$SCRATCH/flye/assembly_info.txt"   "$OUTDIR/" 2>/dev/null || true
cp "$SCRATCH/flye/assembly.fasta"      "$OUTDIR/" 2>/dev/null || true
cp "$SCRATCH/badread.log"              "$OUTDIR/" 2>/dev/null || true
cp "$SCRATCH/flye.stderr"              "$OUTDIR/" 2>/dev/null || true

cat > "$OUTDIR/manifest.json" <<EOF
{
  "isolate": "$ISO", "tag": "$TAG", "arm": "$ARM",
  "depth_fraction_x": $DEPTH, "seed": $SEED,
  "badread_error_model": "$EMODEL", "badread_qscore_model": "$QMODEL",
  "small_plasmid_bias": true,
  "flye_mode": "$FLYEMODE", "threads": $THREADS,
  "read_bases": ${READBASES:-0},
  "seconds_reads": $((t1-t0)), "seconds_flye": $((t2-t1)),
  "flye_version": "$($RUN flye --version 2>&1 | head -1)",
  "badread_version": "$($RUN badread --version 2>&1 | head -1)"
}
EOF

# Reads deleted after graph extraction (proposal Sec. II).
rm -rf "$SCRATCH"
echo "[$TAG] COMPLETE -> $OUTDIR"
