#!/bin/bash

for d in `find *BP_Region* -type d`; 
do 
    echo "Re-plotting $d"
    cd $d
    cp results.txt results.txt.old
    python ../replot.py
    cd ..
done