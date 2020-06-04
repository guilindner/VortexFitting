#!/usr/bin/env/ python
"""
Convert NetCDF4 files to ASCII (plain text)
"""

import argparse
import numpy as np
from netCDF4 import Dataset

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='convert file from netCDF to ASCII format',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', dest='infile',
                        help='input NetCDF file', metavar='FILE')

    parser.add_argument('-o', '--output', dest='outfile',
                        help='output ASCII file', metavar='FILE')

    args = parser.parse_args()

grp1 = Dataset(args.infile, 'r')

u = np.array(grp1.variables['velocity_x'][:, :, :])
v = np.array(grp1.variables['velocity_y'][:, :, :])
w = np.array(grp1.variables['velocity_z'][:, :, :])

print("Converting {:s} file to {:s} file".format(args.infile,args.outfile))

for k in range(len(u)):
    outfile = open(args.outfile, 'w')
    outfile.write('x y u v \n')
    for i in range(u[0, 0].size):
        for j in range(v[0, 0].size):
            outfile.write(str(i)+' '+str(j)+' '+str(u[k, j, i])+' '+str(v[k, j, i])+'\n')

    ###this routine reads the ascii file in top of the netCDF file
    ###used to see if the exported ascii is equal to the original file
    ###put this in vortexfitting.py after "a = VelocityField(args.infilename, args.timestep)"

    #a.uu = []
    #a.vv = []

    #infile = open('ascii/DNS_zPlane0.dat', 'r')
    #lines = infile.readlines()[1:]
    #for x in lines:
    #    a.uu.append(float(x.split(' ')[2]))
    #    a.vv.append(float(x.split(' ')[3]))

    #a.uu = np.array(a.uu)
    #a.vv = np.array(a.vv)
    #a.uu = a.uu.reshape(a.u[:, 0].size, a.u[0, :].size)
    #a.vv = a.vv.reshape(a.v[:, 0].size, a.v[0, :].size)
    #a.u = a.uu
    #a.v = a.vv

