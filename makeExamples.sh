#!/bin/bash

cd `pwd`/src

python vortexfitting.py -i ../data/adim_vel_{:06d}.dat -o ../results/adim_vel_000010 -first 10 -t 5 -b 20 -ft piv_tecplot 

#python vortexfitting.py -i ../data/test_dataHIT.nc -o ../results/test_dataHIT -ft dns

python3 vortexfitting.py -i ../data/test_dataPIV.nc -o ../results/test_dataPIV -ft piv_netcdf -t 1.5

python vortexfitting.py -i ../data/dim_vel_{:06d}.dat -o ../results/dim_vel_000010 -first 10 -mf ../data/mean.dat -ft piv_tecplot -t 50 -rmax 2

#python vortexfitting.py -i ../data/adim_vel_{:06d}.dat -o ../results/adim_vel_000010 -first 10 -last 10 -t 5 -ft piv_tecplot 

