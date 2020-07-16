import os
import sys

is_windows = sys.platform.startswith('win')

if is_windows:
    print("Running on Windows ...")
else:
    print("Running on Linux ...")

cwd = os.getcwd()
os.chdir(cwd + '/tests')

cmd = 'python test_fitting.py'
os.system(cmd)
cmd = 'python test_tools.py'
os.system(cmd)
cmd = 'python testOseen.py'
os.system(cmd)

os.chdir(cwd + '/vortexfitting')

cmd = 'python convertToASCII.py -i ../data/example_data_HIT.nc -o ../data/test_dataHIT_ascii.dat'
os.system(cmd)
cmd = 'python convertToNC.py -i ../data/test_dataHIT_ascii.dat -o ../data/test_dataHIT_back_converted.nc'
os.system(cmd)
cmd = 'python generateNetCDF.py'
os.system(cmd)
