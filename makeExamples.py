import os
import sys

is_windows = sys.platform.startswith('win')

if is_windows:
    print("Running on Windows ...")
else:
    print("Running on Linux ...")

cmd = 'python run.py -i data/example_Ub_planeZ_0.01.raw ' \
      '-ft openfoam -o results/example_openfoam -rmax 0'
os.system(cmd)

cmd = 'python run.py -i data/example_Ub_planeZ_0.01.raw ' \
      '-ft openfoam -xy 20 50 -p detect'
os.system(cmd)

cmd = 'python run.py -i data/example_adim_vel_{:06d}.dat ' \
      '-ft piv_tecplot -o results/example_adim_vel_000010 -first 10 -t 5 -b 20'
os.system(cmd)

cmd = 'python run.py -i data/example_data_HIT.nc ' \
      '-ft dns -o results/example_data_HIT '
os.system(cmd)

cmd = 'python run.py -i data/example_data_numerical_PIV.nc ' \
      '-ft piv_netcdf -o results/example_data_numerical_PIV -t 1.5 -rmax 0'
os.system(cmd)

cmd = 'python run.py -i data/example_dim_vel_{:06d}.dat ' \
      '-ft piv_tecplot -o results/example_dim_vel_000010 -first 10 -mf data/example_mean.dat -t 50 -ct 0.5 -b 15'
os.system(cmd)

cmd = 'python run.py -i data/example_vel_{:06d}.dat ' \
      '-ft piv_tecplot -o results/example_temporal_series -mf data/example_mean.dat -first 5 -last 6 -t 20 -ct 0.5'
os.system(cmd)
