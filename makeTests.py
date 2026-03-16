import os
import sys

PYTHON = sys.executable

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

    cmd = f'{PYTHON} test_fitting.py'
    os.system(cmd)
    cmd = f'{PYTHON} test_tools.py'
    os.system(cmd)
    cmd = f'{PYTHON} test_Lamb_Oseen.py'
    os.system(cmd)
    cmd = f'{PYTHON} test_Rankine.py'
    os.system(cmd)
    
    os.chdir(cwd + '/vortexfitting')

    cmd = f'{PYTHON} convertToASCII.py -i ../data/example_data_HIT.nc -o ../data/test_dataHIT_ascii.dat'
    os.system(cmd)
    cmd = f'{PYTHON} convertToNC.py -i ../data/test_dataHIT_ascii.dat -o ../data/test_dataHIT_back_converted.nc'
    os.system(cmd)

