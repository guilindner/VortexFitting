#!/bin/bash

dir=`pwd`

cd $dir/tests

python test_fitting.py
#python test_schemes.py # TODO test schemes !
python test_tools.py

cd $dir/src

python convertToASCII.py -i ../data/example_dataHIT.nc -o ../data/test_dataHIT_ascii.dat
python convertToNC.py -i ../data/test_dataHIT_ascii.dat -o ../data/test_dataHIT2.nc
python generateNetCDF.py
