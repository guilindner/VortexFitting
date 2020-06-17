#!/bin/bash

dir=`pwd`

cd $dir/tests

python3 test_fitting.py
python3 test_tools.py
python3 testOseen.py
#python3 test_schemes.py # TODO test schemes !

cd $dir/vortexfitting

python3 convertToASCII.py -i ../data/example_dataHIT.nc -o ../data/test_dataHIT_ascii.dat
#OK python3 convertToNC.py -i ../data/test_dataHIT_ascii.dat -o ../data/test_dataHIT_back_converted.nc
#OK python3 generateNetCDF.py
