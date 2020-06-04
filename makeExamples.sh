#!/bin/bash

cd `pwd`/src

#python3 vortexfitting.py -i ../data/example_adim_vel_{:06d}.dat -o ../results/example_adim_vel_000010 -first 10 -t 5 -b 20 -ft piv_tecplot 

#python3 vortexfitting.py -i ../data/example_dataHIT.nc -o ../results/example_dataHIT -ft dns

#python3 vortexfitting.py -i ../data/example_dataPIV.nc -o ../results/example_dataPIV -ft piv_netcdf -t 1.5

python3 vortexfitting.py -i ../data/example_dim_vel_{:06d}.dat -o ../results/example_dim_vel_000010 -first 10 -mf ../data/example_mean.dat -ft piv_tecplot -t 50 -rmax 2

python3 vortexfitting.py -i ../data/example_adim_vel_{:06d}.dat -o ../results/example_adim_vel_000010 -first 10 -last 10 -t 5 -ft piv_tecplot 

