#!/bin/bash
# Put the plasmid env directly on PATH instead of going through
# `micromamba run`. The wrapper forwards signals to its whole process group, so
# a construct like `$(micromamba run flye --version | head -1)` -- where head
# closes the pipe early -- can SIGTERM the calling script. That killed the
# assembly driver twice, ~90 s after the first assembly completed each time.
# Activating by PATH removes the wrapper from the picture entirely (and is
# faster: no per-call process spawn).
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export PLASMID_ENV="$HOME/micromamba/envs/plasmid"
export PATH="$PLASMID_ENV/bin:$PATH"
