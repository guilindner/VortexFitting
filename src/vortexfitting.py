#!/usr/bin/env/ python3
"""vortex detection tool, by Guilherme Lindner, 2017-04\n
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
    parser = argparse.ArgumentParser(description='Optional app description',
                                     formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-i', '--input', dest='infilename',
                        default='../data/test_data.nc',
                        help='input NetCDF file', metavar='FILE')
                        
    parser.add_argument('-o', '--output', dest='outfilename',
                        help='output NetCDF file', metavar='FILE')
    
    parser.add_argument('-s', '--scheme', dest='scheme', type=int, default=2,
                        help='Scheme for differencing\n'
                             '2 = second order\n'
                             '4 = fourth order')
    
    parser.add_argument('-T', '--time', dest='timestep', type=int,
                        default=0,
                        help='Timestep desired')
                        
    parser.add_argument('-d', '--detect', dest='detect',
                        default='swirling',
                        help='Detection method:\n'
                             'Q = Q criterion\n'
                             'swirling = 2D Swirling Strength')
    
    parser.add_argument('-t', '--threshold', dest='threshold',
                        default=0., type=float,
                        help='Threshold for detection, integer')

    parser.add_argument('-b', '--boxsize', dest='boxsize',
                        default=6, type=int,
                        help='Box size for the detection')
    
    args = parser.parse_args()
    
    start = time.time()
    #---- LOAD DATA ----#
    print("Opening file:",args.infilename)
    print("Time:", args.timestep)
    
    a = VelocityField(args.infilename,args.timestep)
    Umean = np.mean(a.u)
    Vmean = np.mean(a.v)
    u2 = a.u - Umean
    v2 = a.v - Vmean
    totalvel = np.sqrt(a.u**2 + a.v**2) 
    
    
    #---- DIFFERENCE APPROXIMATION ----# 
    lap = time.time()
    if args.scheme == 4:
        a.derivative = schemes.fourth_order_diff(a)
    else:
        a.derivative = schemes.second_order_diff(a)
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
    print('x','y','swirl')
    #for i in range(len(peaks[0])):
    #    print(peaks[0][i],peaks[1][i],peaks[2][i])
    #---- PEAKS DIRECTION OF ROTATION ----#
    dirL, dirR = detection.direction_rotation(vorticity,peaks)
    #---- SAVING OUTPUT FILE ----#
    if args.outfilename == None:
        pass
    else:
        print("saving file",args.outfilename)
    #---- PLOTTING ----#
    plt.subplot()
    plt.title('Vorticity, rotation')
    plt.scatter(dirR[0],dirR[1],s=dirR[2]*10,edgecolor='G',facecolor='none')
    plt.scatter(dirL[0],dirL[1],s=dirL[2]*10,edgecolor='Y',facecolor='none')
    plt.imshow(vorticity, interpolation='nearest', cmap="seismic")
    plt.tight_layout()

    xCenter = 1000 
    yCenter = 150
    dist = 20
    X, Y = np.meshgrid(a.dx[xCenter-dist:xCenter+dist],
                       a.dy[yCenter-dist:yCenter+dist])
    U = a.u[yCenter-dist:yCenter+dist,xCenter-dist:xCenter+dist]-Umean
    V = a.v[yCenter-dist:yCenter+dist,xCenter-dist:xCenter+dist]
    
    plt.figure()
    plt.title('Arrows scale with plot width, not view')
    Q = plt.quiver(X, Y, U, V)
 
    #fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)#, sharex=True, sharey=False)
    #ax1.imshow(a.u, cmap='seismic')
    #ax1.set_title('Velocity U (velocity_s)')
    
    #ax2.imshow(a.v, interpolation='nearest', cmap='seismic')
    #ax2.set_title('Velocity V (velocity_n)')
    
    #ax3.set_title('Total velocity')
    #ax3.imshow(totalvel, interpolation='nearest', cmap='seismic')
    #ax3.scatter(peaks[0], peaks[1],s=peaks[2]*10,edgecolors='c',facecolors='none')
    
    #ax4.imshow(vorticity, interpolation='nearest', cmap="seismic")
    #ax4.set_title('Vorticity, rotation')
    #ax4.scatter(dirR[0],dirR[1],s=dirR[2]*10,edgecolors='y',facecolors='y')
    #ax4.scatter(dirL[0],dirL[1],s=dirL[2]*10,edgecolors='g',facecolors='g')
    #plt.tight_layout()
    
    print(round(time.time() - start,3), 'seconds (Total execution time)')
    plt.show()
