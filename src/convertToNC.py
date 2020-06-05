#!/usr/bin/env/ python
"""
Convert ASCII files to NetCDF4 (plain text)
"""

import netCDF4
import argparse
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='convert file from ASCII to netCDF format',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', dest='infile',
                        help='input ASCII file', metavar='FILE')

    parser.add_argument('-o', '--output', dest='outfile',
                        help='output NetCDF file', metavar='FILE')

    parser.add_argument('-nx', '--nx', dest='ndimx',
                        help='spatial mesh dimension, for the x variable', default=159)

    parser.add_argument('-ny', '--ny', dest='ndimy',
                        help='spatial mesh dimension, for the y variable', default=134)

    args = parser.parse_args()

datafile_write = netCDF4.Dataset(args.outfile, 'w', format='NETCDF4')
datafile_write.description = 'Experiments conducted at Rouen ...'

ndimx = args.ndimx  # spacing
ndimy = args.ndimy  # spacing

# dimensions
datafile_write.createDimension('resolution_x', ndimx)
datafile_write.createDimension('resolution_y', ndimy)
datafile_write.createDimension('resolution_z', 1)

# variables
velocity_x = datafile_write.createVariable('velocity_x', 'f4', ('resolution_z',
                                                                'resolution_y',
                                                                'resolution_x'))
velocity_y = datafile_write.createVariable('velocity_y', 'f4', ('resolution_z',
                                                                'resolution_y',
                                                                'resolution_x'))
velocity_z = datafile_write.createVariable('velocity_z', 'f4', ('resolution_z',
                                                                'resolution_y',
                                                                'resolution_x'))
grid_x = datafile_write.createVariable('grid_x', 'f4', 'resolution_x')
grid_y = datafile_write.createVariable('grid_y', 'f4', 'resolution_y')

# data
# velocity_x[:] = np.random.random((1,ndimy,ndimx))/1
# velocity_y[:] = np.random.random((1,ndimy,ndimx))/1
# velocity_z[:] = np.random.random((1,ndimy,ndimx))/1

# grid
x = np.linspace(0, ndimy, ndimx)
y = np.linspace(0, ndimy, ndimx)

print("Converting {:s} file to {:s} file".format(args.infile, args.outfile))

infile = open(args.infile, 'r')
line = infile.readline()
lines = infile.readlines()
for j in range(ndimy):
    for i in range(ndimx):
        velocity_x[0, j, i] = lines[j * ndimx + i].split()[2]
        velocity_y[0, j, i] = lines[j * ndimx + i].split()[3]
        if j == 0:
            grid_x[i] = lines[i].split()[0]
        if i == 0:
            grid_y[j] = lines[j * ndimx].split()[1]

datafile_write.close()
