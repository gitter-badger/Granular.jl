#!/bin/bash

# Optionally use this script to launch profiling.jl with different Julia 
# optimization levels.  Defaults are --optimize=2, --math-mode=ieee.

declare -a arr=(\
    "--procs 1 --optimize=1 --math-mode=ieee" \
    "--procs 1 --optimize=2 --math-mode=ieee" \
    "--procs 1 --optimize=3 --math-mode=ieee" \
    "--procs 1 --optimize=1 --math-mode=fast" \
    "--procs 1 --optimize=2 --math-mode=fast" \
    "--procs 1 --optimize=3 --math-mode=fast")

for flags in "${arr[@]}"; do
    julia --color=yes $flags profiling.jl
done
