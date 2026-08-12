#!/bin/bash
# Chain the remainder: wait for read generation to drain, clean up any Flye
# work-dir left behind by an interrupted run, then assemble everything.
set -uo pipefail

echo "waiting for Badread to drain..."
while pgrep -x badread >/dev/null 2>&1 || pgrep -f 'badread simulate' >/dev/null 2>&1; do
  sleep 20
done
echo "read generation finished at $(date +%H:%M:%S)"
echo "read sets: $(ls -d "$HOME"/scratch/*/ 2>/dev/null | wc -l)"

# an interrupted Flye leaves a partial work dir that would make it resume wrongly
pkill -f flye-modules 2>/dev/null
sleep 3
for d in "$HOME"/scratch/*/flye; do
  [ -d "$d" ] || continue
  tag=$(basename "$(dirname "$d")")
  if [ ! -s "/mnt/d/DNA-Sequencing/work/asm/$tag/assembly_graph.gfa" ]; then
    echo "  clearing partial flye dir for $tag"
    rm -rf "$d"
  fi
done

echo "starting Phase B at $(date +%H:%M:%S)"
# No JOBS / THREADS override. 08 is sequential by design and pins threads per
# arm (6 r9_raw, 10 r10_hq) because Flye's peak memory scales with thread count
# and the R9 arm OOM-killed this VM three times. The old `JOBS=2 THREADS=11`
# here was dead (08 reads neither) but it read as an instruction to parallelise,
# which is exactly the change that must not be made.
bash /mnt/d/DNA-Sequencing/work/scripts/08_batch_flye.sh
echo "ALL DONE at $(date +%H:%M:%S)"
