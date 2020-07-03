#!/usr/bin/env/ python3
"""
Generate a netCDF file with a vortex field.
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

import fitting

args = []

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='generate a vortex field in a netCDF file',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-o', '--output', dest='outfile', type=str,
                        help='output NetCDF file', metavar='FILE',
                        default='../data/generatedField.nc')

    parser.add_argument('-ndim', '--ndim', dest='ndim', type=int,
                        help='spatial mesh dimension, for each x and y variables',
                        default=256)

    args = parser.parse_args()

print('Generating {:s} file with a {:d}x{:d} mesh'.format(args.outfile, args.ndim, args.ndim))

# Try to write the file
try:
    datafile_write = netCDF4.Dataset(args.outfile, 'w', format='NETCDF4')
except IOError:
    print('There was an error writing the file!')
    sys.exit()

datafile_write.description = 'Sample field with an Oseen vortex'

ndim = args.ndim  # spacing

# dimensions
datafile_write.createDimension('resolution_x', ndim)
datafile_write.createDimension('resolution_y', ndim)
datafile_write.createDimension('resolution_z', 1)

# variables
velocity_x = datafile_write.createVariable('velocity_x', 'f4', ('resolution_z', 'resolution_y', 'resolution_x'))
velocity_y = datafile_write.createVariable('velocity_y', 'f4', ('resolution_z', 'resolution_y', 'resolution_x'))
velocity_z = datafile_write.createVariable('velocity_z', 'f4', ('resolution_z', 'resolution_y', 'resolution_x'))

# data
velocity_x[:] = np.random.random((1, ndim, ndim)) / 10
velocity_y[:] = np.random.random((1, ndim, ndim)) / 10
velocity_z[:] = np.random.random((1, ndim, ndim)) / 10

# grid
x_grid = np.linspace(0, ndim, ndim)
y_grid = np.linspace(0, ndim, ndim)

x_matrix, y_matrix = np.meshgrid(x_grid, y_grid)
core_radius = 5.0
gamma = 30
x_center = 64
y_center = 192
u_advection = 0.0
v_advection = 0.3
u_data, v_data = fitting.velocity_model(core_radius, gamma, x_center, y_center, u_advection, v_advection, x_matrix,
                                        y_matrix)
u_data = u_data + u_advection
v_data = v_data + v_advection

velocity_x[0, :, :] += u_data[:, :]
velocity_y[0, :, :] += v_data[:, :]
s = 4  # sampling factor for quiver plot
plt.quiver(x_matrix[::s, ::s], y_matrix[::s, ::s], velocity_x[0, ::s, ::s], velocity_y[0, ::s, ::s])

plt.show()
datafile_write.close()
