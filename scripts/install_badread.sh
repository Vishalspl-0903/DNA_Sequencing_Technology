#!/bin/bash
set -e
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"
"$MM" install -y -n plasmid -r "$MAMBA_ROOT_PREFIX" -c conda-forge -c bioconda badread
echo "--- verify all ---"
for t in flye minimap2 seqkit rasusa badread; do
  printf '%-10s ' "$t"
  "$MM" run -n plasmid -r "$MAMBA_ROOT_PREFIX" bash -c "command -v $t >/dev/null 2>&1 && ($t --version 2>&1 | head -1) || echo MISSING"
done
