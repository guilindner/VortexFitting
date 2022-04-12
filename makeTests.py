import os
import sys

is_windows = sys.platform.startswith('win')

if is_windows:
    print("Running on Windows ...")
else:
    print("Running on Linux ...")

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")
else:        
    cwd = os.getcwd()
    os.chdir(cwd + '/tests')    
    cwd = os.getcwd()
    os.chdir(cwd + '/tests')

    cmd = 'python3 test_fitting.py'
    os.system(cmd)
    cmd = 'python3 test_tools.py'
    os.system(cmd)
    cmd = 'python3 testOseen.py'
    os.system(cmd)

    os.chdir(cwd + '/vortexfitting')

    cmd = 'python3 convertToASCII.py -i ../data/example_data_HIT.nc -o ../data/test_dataHIT_ascii.dat'
    os.system(cmd)
    cmd = 'python3 convertToNC.py -i ../data/test_dataHIT_ascii.dat -o ../data/test_dataHIT_back_converted.nc'
    os.system(cmd)
    cmd = 'python3 generateNetCDF.py'
    os.system(cmd)
