#!/bin/bash
# minimap2 lives in the plasmid env, so the labeller runs under WSL.
# Env activated by PATH, not `micromamba run` (see env.sh for why).
set -uo pipefail
source /mnt/d/DNA-Sequencing/work/scripts/env.sh
export PIPE_ROOT=/mnt/d/DNA-Sequencing/work
python /mnt/d/DNA-Sequencing/work/scripts/04_label.py "$(command -v minimap2)"
