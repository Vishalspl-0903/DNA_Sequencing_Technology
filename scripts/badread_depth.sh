#!/bin/bash
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"
SRC=$("$MM" run -n plasmid -r "$MAMBA_ROOT_PREFIX" python -c "import badread,os;print(os.path.dirname(badread.__file__))")
echo "source: $SRC"
echo "=========== depth handling ==========="
grep -rn "depth" "$SRC"/*.py | head -40
echo "=========== small_plasmid_bias ==========="
grep -rn "small_plasmid_bias" "$SRC"/*.py | head -20
