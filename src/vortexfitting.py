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
import identification
from classes import VelocityField

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Optional app description',formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-i', '--input', dest='infilename', default='../data/test_data.nc',
                        help='input NetCDF file', metavar='FILE')
                        
    parser.add_argument('-o', '--output', dest='outfilename',
                        help='input NetCDF file', metavar='FILE')
    
    parser.add_argument('-s', '--scheme', dest='scheme', type=int, default=2,
                        help='Scheme for differencing:\n2 = second order\n4 = fourth order')
    
    parser.add_argument('-T', '--time', dest='timestep', type=int, default=0,
                        help='Timestep desired')
                        
    parser.add_argument('-d', '--detect', dest='detect', default='swirling',
                        help='Detection method:\nQ = Q criterion\nswirling = 2D Swirling Strength')
    
    parser.add_argument('-t', '--threshold', dest='threshold', type=float, default=0.,
                        help='Threshold for detection, integer')

    parser.add_argument('-b', '--boxsize', dest='boxsize', type=int, default=6,
                        help='Box size for the detection')
    
    args = parser.parse_args()
    
    start = time.time()
    #---- LOAD DATA ----#
    print("Opening file:",args.infilename)
    print("Time:", args.timestep)
    
    a = VelocityField(args.infilename,args.timestep)
    totalvel = np.sqrt(a.u**2 + a.v**2) 
    
    
    #---- DIFFERENCE APPROXIMATION ----# 
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
    
    #---- METHOD FOR DETECTION OF VORTICES ----#
    lap = time.time()
    if args.detect == 'Q':
        detected = identification.q_criterion(a)
    elif args.detect == 'swirling':
        detected = identification.calc_swirling(a)
    print(round(time.time() - lap,3), 'seconds')

    #---- PEAK DETECTION ----#
    print("Detecting peak of swirling strength")
    print("threshold=",args.threshold,"box size=",args.boxsize)

    peaks = detection.find_peaks(detected, args.threshold, args.boxsize)
    print("Vortices found:",len(peaks[0]))

    #---- PEAKS DIRECTION OF ROTATION ----#
    clockwise, counterclockwise = detection.direction_rotation(vorticity,peaks)

    #---- SAVING OUTPUT FILE ----#
    if args.outfilename == None:
        pass
    else:
        print("saving file",args.outfilename)
    
    #---- PLOTTING ----#
    
    #~ X, Y = np.meshgrid(a.dx[0:15],a.dy[180:230])
    #~ U = a.u[0:15,180:230]
    #~ V = a.v[0:15,180:230]
    #~ #print(U)
    #~ plt.figure()
    #~ plt.title('Arrows scale with plot width, not view')
    #~ Q = plt.quiver(X, Y, U, V,pivot='mid')
 
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)#, sharex=True, sharey=False)
    ax1.imshow(a.u, cmap='seismic')
    ax1.set_title('Velocity U (velocity_s)')
    
    ax2.imshow(a.v, interpolation='nearest', cmap='seismic')
    ax2.set_title('Velocity V (velocity_n)')
    
    ax3.set_title('Total velocity')
    #ax3.set_xlim(0,1152)
    #ax3.set_ylim(250,0)
    ax3.imshow(totalvel, interpolation='nearest', cmap='seismic')
    ax3.scatter(peaks[0], peaks[1],s=peaks[2],edgecolors='c',facecolors='none')
    
    ax4.imshow(vorticity, interpolation='nearest', cmap="seismic")
    ax4.set_title('Vorticity, rotation')
    #ax4.set_xlim(0,1152)
    #ax4.set_ylim(250,0)
    ax4.scatter(counterclockwise[0],counterclockwise[1],counterclockwise[2],facecolors='y')
    ax4.scatter(clockwise[0],clockwise[1],clockwise[2],facecolors='g')
    plt.tight_layout()
    
    print(round(time.time() - start,3), 'seconds (Total execution time)')
    plt.show()
    
    
    
    
