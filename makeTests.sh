#!/bin/bash

dir=`pwd`

cd $dir/tests

python3 test_fitting.py
python3 test_tools.py
python3 testOseen.py

cd $dir/vortexfitting

python3 convertToASCII.py -i ../data/example_data_HIT.nc -o ../data/test_dataHIT_ascii.dat
python3 convertToNC.py -i ../data/test_dataHIT_ascii.dat -o ../data/test_dataHIT_back_converted.nc
python3 generateNetCDF.py
