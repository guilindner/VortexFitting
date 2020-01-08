#!/bin/bash
python3 vortexfitting.py -i ../data/f026_vr005_vs20/Ub_PlaneZ_35mm_{:1d}.raw -mf ../data/f026_vr005_vs20/UbMean_planeZ_35mm.raw -t 10 -rmax 0.001 -first 1 -last 30 -step 1 -o ../results/f026_vr005_vs20/
convert -delay 20 -loop 0 ../results/f026_vr005_vs20/accepted_{1..30..1}.svg ../results/f026_vr005_vs20/animated.gif
