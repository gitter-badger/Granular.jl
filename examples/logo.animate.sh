#!/bin/bash
convert -trim +repage \
    -delay 10 -transparent-color white -reverse -loop 0 \
    logo/*.png logo/logo.gif
