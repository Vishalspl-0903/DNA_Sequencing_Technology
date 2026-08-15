#!/bin/bash
# Userspace bioinformatics toolchain via micromamba (no sudo required).
# Note: this Ubuntu image has no bzip2, so the micromamba tarball is
# extracted with python3's tarfile module instead of tar.
set -e
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"

mkdir -p "$HOME/bin"

if [ ! -x "$MM" ]; then
  echo "[1/4] downloading micromamba..."
  curl -Ls -o /tmp/micromamba.tar.bz2 \
    https://micro.mamba.pm/api/micromamba/linux-64/latest
  ls -l /tmp/micromamba.tar.bz2
  python3 - <<'PY'
import tarfile, shutil, os
with tarfile.open("/tmp/micromamba.tar.bz2", "r:bz2") as t:
    m = t.extractfile("bin/micromamba")
    dest = os.path.expanduser("~/bin/micromamba")
    with open(dest, "wb") as o:
        shutil.copyfileobj(m, o)
    os.chmod(dest, 0o755)
print("extracted ->", dest)
PY
fi
echo "micromamba: $($MM --version)"

echo "[2/4] creating env 'plasmid' (flye, minimap2, seqkit, rasusa, python 3.11)..."
"$MM" create -y -n plasmid -r "$MAMBA_ROOT_PREFIX" \
  -c conda-forge -c bioconda \
  python=3.11 flye minimap2 seqkit rasusa pip

echo "[3/4] installing Badread (pip, into env)..."
"$MM" run -n plasmid -r "$MAMBA_ROOT_PREFIX" pip install --quiet badread

echo "[4/4] verifying..."
for t in flye minimap2 seqkit rasusa badread; do
  printf '%-10s ' "$t"
  "$MM" run -n plasmid -r "$MAMBA_ROOT_PREFIX" bash -c "command -v $t >/dev/null 2>&1 && ($t --version 2>&1 | head -1) || echo MISSING"
done
echo "DONE"
