#!/bin/bash

# Status script which traverses the subdirectories of the current folder for 
# simulations.  You may want to add this to your shell's PATH variable.

set -e
cmd='julia --color=yes -e "import SeaIce; SeaIce.status()"'

if [[ "$1" == "loop" ]]; then
    while true; do
        date
        eval $cmd
        sleep 10
    done
else
    eval $cmd
fi
