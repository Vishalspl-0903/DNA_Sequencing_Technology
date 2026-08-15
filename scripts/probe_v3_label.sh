#!/bin/bash
# Label the probe assemblies with the SAME unit-tested labeller that produced
# the cohort labels, pointed at work/probe via PIPE_ROOT. Reusing 04_label.py
# rather than reimplementing it is the point: the probe's answer has to be
# commensurable with the cohort's.
set -uo pipefail
source /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/scripts/env.sh
export PIPE_ROOT=/mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/probe
python /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/scripts/04_label.py \
       "$(command -v minimap2)"
