#!/bin/bash
# Label the round-3 probe with the same unit-tested labeller, pointed at
# work/probe3. Separate file rather than an argument because wsl.exe does not
# reliably preserve a quoted `bash -c '...'` string across the interop boundary.
set -uo pipefail
source /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/scripts/env.sh
export PIPE_ROOT=/mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/probe3
python /mnt/d/Plasmid-GNN/DNA_Sequencing_Technology/work/scripts/04_label.py \
       "$(command -v minimap2)"
