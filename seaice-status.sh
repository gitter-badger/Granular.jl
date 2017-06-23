#!/bin/bash

# Status script which traverses the subdirectories of the current folder for 
# simulations.  You may want to add this to your shell's PATH variable.

set -e
cmd_sing='julia --color=yes -e "import SeaIce; SeaIce.status()"'
cmd_loop='julia --color=yes -e "import SeaIce; SeaIce.status(loop=true, t_int=10)"'

if [[ "$1" == "loop" ]]; then
    eval $cmd_loop
else
    eval $cmd_sing
fi
