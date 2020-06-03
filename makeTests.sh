#!/bin/bash

cd `pwd`/tests

python test_fitting.py
#python test_schemes.py
python test_tools.py

cd `pwd`/src

python convertToASCII.py -i ../data/test_dataHIT.nc -o ../data/test_dataHIT_ascii.dat
python convertToNC.py -i ../data/test_dataHIT_ascii.dat -o ../data/test_dataHIT2.nc
python generateNetCDF.py
