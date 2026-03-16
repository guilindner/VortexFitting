#!/usr/bin/env python3
"""
Convert NetCDF4 files to ASCII (plain text)
"""

import sys
import argparse
import numpy as np
import netCDF4


def netcdf_to_ascii(input_file: str, output_file: str) -> None:
    """
    Convert a NetCDF file to an ASCII file.

    :param input_file: path to the input NetCDF file
    :type input_file: str
    :param output_file: path to the output ASCII file
    :type input_file: str

    :returns: empty
    :rtype: empy
    """
    # Try to read the NetCDF file
    try:
        datafile_read = netCDF4.Dataset(input_file, 'r')
    except IOError:
        print(f"There was an error opening the file: {input_file}")
        sys.exit()

    # Extract velocity data
    try:
        u = np.array(datafile_read.variables['velocity_x'][:, :, :])
        v = np.array(datafile_read.variables['velocity_y'][:, :, :])
        w = np.array(datafile_read.variables['velocity_z'][:, :, :])
    except KeyError as e:
        print(f"Missing expected variable in NetCDF file: {e}")
        sys.exit()

    print(f"Converting {input_file} file to {output_file} file")
    # Try to write to the ASCII file
    try:
        with open(output_file, 'w') as outfile:
            for k in range(len(u)):
                outfile.write('x y u v \n')
                for i in range(u.shape[2]):
                    for j in range(u.shape[1]):
                        outfile.write(f"{i} {j} {u[k, j, i]} {v[k, j, i]}\n")
    except IOError:
        print(f"There was an error writing the file: {output_file}")
        sys.exit()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert file from NetCDF to ASCII format',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', dest='infile', type=str,
                        default='../data/test_netCDF.nc',
                        help='Input NetCDF file', metavar='FILE')

    parser.add_argument('-o', '--output', dest='outfile', type=str,
                        default='../data/test_netCDF_ascii.dat',
                        help='Output ASCII file', metavar='FILE')

    args = parser.parse_args()

    # Call the function with the provided arguments
    netcdf_to_ascii(args.infile, args.outfile)

