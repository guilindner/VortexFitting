#!/bin/bash

cd `pwd`/src

python vortexfitting.py -i ../data/test_dataHIT.nc -o ../results/test_dataHIT -ft dns

python vortexfitting.py -i ../data/test_dataPIV.nc -o ../results/test_dataPIV -ft piv

python vortexfitting.py -i ../data/dim_vel_{:06d}.dat -o ../results/dim_vel_000010 -first 10 -last 10 -mf ../data/mean.dat -ft tecplot -t 50 -rmax 1

python vortexfitting.py -i ../data/adim_vel_{:06d}.dat -o ../results/adim_vel_000010 -first 10 -last 10 -ft tecplot -t 5

