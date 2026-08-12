#!/bin/bash
# Robust, restartable assembly driver.
#
# Phase B died once with the WSL VM going down while two Flye jobs were in
# flight, so this runs SEQUENTIALLY (JOBS=1) with more threads per job instead.
# Everything here is idempotent: a tag whose GFA already exists is skipped, and
# a partial Flye work dir left by an interrupted run is cleared rather than
# resumed (Flye would otherwise try to continue from an inconsistent state).
#
# The outer loop retries so a single termination does not lose the batch; it
# stops when a pass makes no progress.
set -uo pipefail

DROOT="/mnt/d/DNA-Sequencing/work"
MAXPASS="${MAXPASS:-6}"

count_graphs () { ls -d "$DROOT"/asm/*/assembly_graph.gfa 2>/dev/null | wc -l; }
count_reads  () { ls -d "$HOME"/scratch/*/ 2>/dev/null | wc -l; }

for pass in $(seq 1 "$MAXPASS"); do
  before=$(count_graphs)
  remaining=$(count_reads)
  echo "=== pass $pass : graphs=$before  read_sets_left=$remaining  $(date +%H:%M:%S) ==="
  if [ "$remaining" -eq 0 ]; then
    echo "no read sets left to assemble"
    break
  fi

  # clear partial work dirs from any interrupted run
  for d in "$HOME"/scratch/*/flye; do
    [ -d "$d" ] || continue
    tag=$(basename "$(dirname "$d")")
    if [ ! -s "$DROOT/asm/$tag/assembly_graph.gfa" ]; then
      echo "  clearing partial flye dir: $tag"
      rm -rf "$d"
    fi
  done

  # sequential; threads set per arm inside 08 (r9_raw is memory-hungry)
  bash "$DROOT/scripts/08_batch_flye.sh"

  after=$(count_graphs)
  echo "=== pass $pass done : graphs $before -> $after ==="
  if [ "$after" -le "$before" ]; then
    echo "no progress in pass $pass, stopping"
    break
  fi
done

echo "ASSEMBLY DRIVER FINISHED graphs=$(count_graphs) read_sets_left=$(count_reads) $(date +%H:%M:%S)"
