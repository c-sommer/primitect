#!/bin/bash

while read -r -a line
do
    stamp=${line[0]}
    depth_png="_depth.png"
    rgb_jpg="_rgb.jpg"
    
    unzip $stamp.zip -d $stamp
    cp $stamp/depth/0000001-000000000000.png $stamp$depth_png
    cp $stamp/rgb/0000001-000000000000.jpg $stamp$rgb_jpg
    
    rm -r -f $stamp
    rm $stamp.zip
    echo $stamp
done < stamps.dat

echo "Done."
