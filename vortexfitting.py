#!/usr/bin/env/ python
"""vortex detection tool, by Guilherme Lindner, 2017-04
This program load NetCDF files from DNS simulations  or PIV experiments
and detect the vortices and apply a fitting to them.
"""
import sys
import argparse
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy import ndimage
import detection
import schemes
from classes import VelocityField


parser = argparse.ArgumentParser(description='Optional app description',formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-i', '--input', dest='infilename', default='Data/test_data.nc',
                    help='input NetCDF file', metavar='FILE')
                    
parser.add_argument('-o', '--output', dest='outfilename',
                    help='input NetCDF file', metavar='FILE')

parser.add_argument('-s', '--scheme', dest='scheme', type=int, default=2,
                    help='Scheme for differencing:\n2 = second order\n4 = fourth order')

parser.add_argument('-t', '--time', dest='timestep', type=int, default=0,
                    help='Timestep desired')
                    
parser.add_argument('-d', '--detect', dest='detect', type=int, default=0,
                    help='Detection method:\n2 = 2D Swirling Strength')

args = parser.parse_args()

#---- LOAD DATA ----#
print("Opening file:",args.infilename)
print("Time:", args.timestep)

a = VelocityField(args.infilename,args.timestep)
 

#print('Resolution x,y:',a.sizex,a.sizey)

#---- Second-order difference approximation ----#

if args.scheme == 4:
    a.derivative = schemes.fourthOrderDiff(a)
else:
    a.derivative = schemes.secondOrderDiff(a)
strength = []


#---- VORTICITY ----#
print("Calculating vorticity")
vorticity = a.derivative['dudy'] - a.derivative['dvdx']

#---- VELOCITY TENSOR GRADIENT ----#
print("Calculating Velocity Tensor Gradient")
L = np.zeros((a.sizex,a.sizey))
ls2 = np.zeros((a.sizex,a.sizey))

for i in range(a.sizex):
        for j in range(a.sizey):
            A = [[a.derivative['dudx'][i,j],a.derivative['dudy'][i,j],
            a.derivative['dudz'][i,j]],[a.derivative['dvdx'][i,j],
            a.derivative['dvdy'][i,j],a.derivative['dvdz'][i,j]],
            [a.derivative['dwdx'][i,j],a.derivative['dwdy'][i,j],
            -a.derivative['dudx'][i,j]-a.derivative['dvdy'][i,j]]]
            ls = np.linalg.eigvals(A)
            ls2[i,j] = max(abs(ls.imag))


maxima1 = ndimage.maximum_filter(ls2,size=(6,6))
maxima = detection.find_peaks(ls2, 1.0, box_size=6)
strength = maxima[2]*10

if args.outfilename == None:
    pass
else:
    print("saving file",args.outfilename)
    


#---- PLOTTING ----#
# Make plot with vertical (default) colorbar
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
cax = ax1.imshow(a.u, interpolation='nearest', cmap=plt.cm.coolwarm)
ax1.set_title('Velocity U (velocity_s)')
#cbar = fig.colorbar(cax)#, ticks=[-1, 0, 1])

cax = ax2.imshow(a.v, interpolation='nearest', cmap=plt.cm.coolwarm)
ax2.set_title('Velocity V (velocity_n)')
#cbar = fig.colorbar(cax)#, ticks=[-1, 0, 1])

#cax = ax3.imshow(vorticity, interpolation='nearest', cmap=plt.cm.coolwarm)
ax3.set_title('Swirling Strength, max filter')
ax3.set_xlim(0,1152)
ax3.imshow(maxima1)
ax3.scatter(maxima[0], maxima[1],s=strength,edgecolors='r',facecolors='none')

cax = ax4.imshow(ls2, interpolation='nearest', cmap="Greys")
ax4.set_title('Swirling Strength, no filter')
#ax4.scatter(localy,localx,s=strength,edgecolors='r',facecolors='none')
#ax4.plot(maxima['x_peak'], maxima['y_peak'], ls='none',marker='o')
ax4.scatter(maxima[0], maxima[1],s=maxima[2],edgecolors='r',facecolors='none')
#ax4.scatter(maxima['x_peak'], maxima['y_peak'],s=maxima['peak_value'],edgecolors='r',facecolors='none')

plt.tight_layout()
plt.show()
