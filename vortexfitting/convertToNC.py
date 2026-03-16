#!/usr/bin/env python3
"""
Convert ASCII files to NetCDF4 (plain text)
"""

import sys
import argparse
# import numpy as np
import netCDF4


def ascii_to_netcdf(input_file: str, output_file: str) -> None:
    """
    Convert an ASCII file to a NetCDF file.

s    :param input_file: path to the input ASCII file
    :type input_file: str
    :param output_file: path to the output NetCDF file
    :type input_file: str    
    :returns: empty
    :rtype: empy    
    """
    ndimx, ndimy = deduce_dimensions(input_file)
    # Try to create the NetCDF file
    try:
        datafile_write = netCDF4.Dataset(output_file, 'w', format='NETCDF4')
    except IOError:
        print(f"There was an error writing the file: {output_file}")
        sys.exit()

    datafile_write.description = 'Experiments conducted at Rouen ...'

    # Define dimensions
    datafile_write.createDimension('resolution_x', ndimx)
    datafile_write.createDimension('resolution_y', ndimy)
    datafile_write.createDimension('resolution_z', 1)

    # Define variables
    velocity_x = datafile_write.createVariable('velocity_x', 'f4', ('resolution_z', 'resolution_y', 'resolution_x'))
    velocity_y = datafile_write.createVariable('velocity_y', 'f4', ('resolution_z', 'resolution_y', 'resolution_x'))
    # Uncomment for velocity_z
    # velocity_z = datafile_write.createVariable('velocity_z', 'f4', ('resolution_z', 'resolution_y', 'resolution_x'))
    grid_x = datafile_write.createVariable('grid_x', 'f4', 'resolution_x')
    grid_y = datafile_write.createVariable('grid_y', 'f4', 'resolution_y')

    # Grid data
    # x = np.linspace(0, ndimy, ndimx)
    # y = np.linspace(0, ndimy, ndimx)

    print(f"Converting {input_file} file to {output_file} file")

    # Try to read the ASCII file
    try:
        with open(input_file, 'r') as infile:
            infile.readline()  # Skip the header
            lines = infile.readlines()
    except IOError:
        print(f"There was an error reading the file: {input_file}")
        sys.exit()

    # Parse data and populate NetCDF variables
    for j in range(ndimy):
        for i in range(ndimx):
            line_index = j * ndimx + i
            values = lines[line_index].split()
            velocity_x[0, j, i] = float(values[2])
            velocity_y[0, j, i] = float(values[3])
            if j == 0:
                grid_x[i] = float(values[0])
            if i == 0:
                grid_y[j] = float(values[1])

    # Close the NetCDF file
    datafile_write.close()



def deduce_dimensions(input_file: str) -> tuple[int, int]:
    """
    Deduce the dimensions ndimx and ndimy from the input ASCII file.

    :param input_file: path to the input ASCII file
    :type input_file: str
    :returns: ndimx, ndimy, number of points in the x /ydimension
    :rtype: int, int
    """
    with open(input_file, 'r') as infile:
        infile.readline()  # Skip the header line
        first_line = infile.readline()
        # x_first = float(first_line.split()[0])  # First x value
        # y_first = float(first_line.split()[1])  # First y value

        # Read the rest of the file
        lines = infile.readlines()

        # Count unique x and y values
        x_values = {float(line.split()[0]) for line in lines}
        y_values = {float(line.split()[1]) for line in lines}

        ndimx = len(x_values)
        ndimy = len(y_values)

        return ndimx, ndimy


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert file from ASCII to NetCDF format',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', dest='infile', type=str,
                        default='../data/test_dataHIT_ascii.dat',
                        help='Input ASCII file', metavar='FILE')

    parser.add_argument('-o', '--output', dest='outfile', type=str,
                        default='../data/test_dataHIT_back_converted.nc',
                        help='Output NetCDF file', metavar='FILE')

    parser.add_argument('-nx', '--nx', dest='ndimx', type=int,
                        help='Spatial mesh dimension, for the x variable', default=159)

    parser.add_argument('-ny', '--ny', dest='ndimy', type=int,
                        help='Spatial mesh dimension, for the y variable', default=134)

    args = parser.parse_args()

    # Call the function with the provided arguments
    ascii_to_netcdf(args.infile, args.outfile)

