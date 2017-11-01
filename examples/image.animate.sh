#!/bin/bash
convert -trim +repage \
    -delay 10 -transparent-color white -reverse -loop 0 \
    image/*.png image/image.gif

convert -trim +repage \
    -delay 10 -transparent-color white -loop 0 \
    image/*.png image/image-forward.gif
