#!/bin/bash
# minimap2 lives in the plasmid env, so the labeller runs under WSL.
# Env activated by PATH, not `micromamba run` (see env.sh for why).
set -uo pipefail
source /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/scripts/env.sh
export PIPE_ROOT=/mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work
python /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/scripts/04_label.py "$(command -v minimap2)"
