#!/usr/bin/env/ python
"""vortex detection tool, by Guilherme Lindner, 2017-04
This program load NetCDF files from DNS simulations  or PIV experiments
and detect the vortices and apply a fitting to them.
"""
import sys
import argparse
import time
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy import ndimage
import detection
import schemes
from classes import VelocityField

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Optional app description',formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-i', '--input', dest='infilename', default='../data/test_data.nc',
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
    
    start = time.time()
    
    #---- LOAD DATA ----#
    print("Opening file:",args.infilename)
    print("Time:", args.timestep)
    
    a = VelocityField(args.infilename,args.timestep)
     
    
    #print('Resolution x,y:',a.sizex,a.sizey)
    
    #---- Second-order difference approximation ----#
    lap = time.time()
    if args.scheme == 4:
        a.derivative = schemes.fourthOrderDiff(a)
    else:
        a.derivative = schemes.secondOrderDiff(a)
    print(round(time.time() - lap,3), 'seconds') 

    strength = []
    
    #---- VORTICITY ----#
    print("Calculating vorticity")
    vorticity = a.derivative['dudy'] - a.derivative['dvdx']
    
    #---- VELOCITY TENSOR GRADIENT ----#
    print("Calculating Velocity Tensor Gradient")
    lap = time.time()
    A = np.zeros((a.sizex*a.sizey,3,3))

    A = np.array([[a.derivative['dudx'].ravel(),a.derivative['dudy'].ravel(),
                a.derivative['dudz'].ravel()],[a.derivative['dvdx'].ravel(),
                a.derivative['dvdy'].ravel(),a.derivative['dvdz'].ravel()],
                [a.derivative['dwdx'].ravel(),a.derivative['dwdy'].ravel(),
                -a.derivative['dudx'].ravel()-a.derivative['dvdy'].ravel()]])
    A = A.transpose(2,1,0)
    eigenvalues = np.linalg.eigvals(A)
    swirling = np.max(eigenvalues.imag,axis=1).reshape(a.sizex,a.sizey)
       
    print(round(time.time() - lap,3), 'seconds')            
    
    maxima1 = ndimage.maximum_filter(swirling,size=(6,6))
    maxima = detection.find_peaks(swirling, 1.0, box_size=6)
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
    
    cax = ax4.imshow(swirling, interpolation='nearest', cmap="Greys")
    ax4.set_title('Swirling Strength, no filter')
    #ax4.scatter(localy,localx,s=strength,edgecolors='r',facecolors='none')
    #ax4.plot(maxima['x_peak'], maxima['y_peak'], ls='none',marker='o')
    ax4.scatter(maxima[0], maxima[1],s=maxima[2],edgecolors='r',facecolors='none')
    #ax4.scatter(maxima['x_peak'], maxima['y_peak'],s=maxima['peak_value'],edgecolors='r',facecolors='none')
    
    plt.tight_layout()
    print(round(time.time() - start,3), 'seconds (Total execution time)')
    plt.show()
    
    
