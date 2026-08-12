#!/bin/bash
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"
MM2=$("$MM" run -n plasmid -r "$MAMBA_ROOT_PREFIX" bash -c 'command -v minimap2')
echo "minimap2: $MM2"
"$MM" run -n plasmid -r "$MAMBA_ROOT_PREFIX" python \
   /mnt/d/DNA-Sequencing/work/scripts/test_labeller.py "$MM2"
