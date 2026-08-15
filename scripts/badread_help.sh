#!/bin/bash
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"
"$MM" run -n plasmid -r "$MAMBA_ROOT_PREFIX" badread simulate --help 2>&1 | head -60
echo "=========== error models available ==========="
"$MM" run -n plasmid -r "$MAMBA_ROOT_PREFIX" python -c "
import badread, os
d = os.path.join(os.path.dirname(badread.__file__), 'error_models')
q = os.path.join(os.path.dirname(badread.__file__), 'qscore_models')
print('error_models:', sorted(os.listdir(d)) if os.path.isdir(d) else 'n/a')
print('qscore_models:', sorted(os.listdir(q)) if os.path.isdir(q) else 'n/a')
"
