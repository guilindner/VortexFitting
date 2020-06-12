#!/bin/bash

python3 vortexfitting/vortexfitting.py -i data/example_Ub_planeZ_0.01.raw -o results/example_openfoam -ft openfoam -rmax 0

python3 vortexfitting/vortexfitting.py -i data/example_Ub_planeZ_0.01.raw -ft openfoam -xy 20 50

python3 vortexfitting/vortexfitting.py -i data/example_adim_vel_{:06d}.dat -o results/example_adim_vel_000010 -first 10 -t 5 -b 20 -ft piv_tecplot 

python3 vortexfitting/vortexfitting.py -i data/example_dataHIT.nc -o results/example_dataHIT -ft dns

python3 vortexfitting/vortexfitting.py -i data/example_dataPIV.nc -o results/example_dataPIV -ft piv_netcdf -t 1.5

python3 vortexfitting/vortexfitting.py -i data/example_dim_vel_{:06d}.dat -o results/example_dim_vel_000010 -first 10 -mf data/example_mean.dat -ft piv_tecplot -t 50 -rmax 2

python3 vortexfitting/vortexfitting.py -i data/example_adim_vel_{:06d}.dat -o results/example_adim_vel_000010 -first 10 -last 10 -t 5 -ft piv_tecplot 

python3 vortexfitting/vortexfitting.py -i data/example_vel_{:06d}.dat -o results/example_temporal_series -mf data/example_mean.dat -first 5 -last 6 -t 50 -ct 0.5 -ft piv_tecplot 

