#!/usr/bin/env/ python3
"""vortex detection tool, by Guilherme Lindner, 2017-04\n
This program load NetCDF files from DNS simulations  or PIV experiments
and detect the vortices and apply a fitting to them.
"""
import sys
import argparse
import time
import numpy as np

import tools
import plot
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
                        help='Timestep/Sample desired')
                        
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
    
    parser.add_argument('-p', '--plot', dest='plot_x',
                        default='detect',
                        help='Plot on screen:\n'
                             'detect = Vortices position\n'
                             'fields = Velocity fields\n'
                             'quiver = Vector on specific position')
    
    args = parser.parse_args()
    
    start = time.time()
    #---- LOAD DATA ----#
    print("Opening file:",args.infilename)

    print("Sample target: (todo)", args.timestep)
    
    a = VelocityField(args.infilename,args.timestep)
    print("Samples:", a.samples)

    #---- DIFFERENCE APPROXIMATION ----# 
    lap = time.time()
    if args.scheme == 4:
        a.derivative = schemes.fourth_order_diff(a)
    elif args.scheme == 2:
        a.derivative = schemes.second_order_diff(a)
    else:
        print('No scheme', args.scheme, 'found. Exitting!')
        sys.exit()
    print(round(time.time() - lap,3), 'seconds') 
    
    #---- VORTICITY ----#
    print("Calculating vorticity")

    vorticity = a.derivative['dvdx'] - a.derivative['dudy']

    #---- METHOD FOR DETECTION OF VORTICES ----#
    lap = time.time()
    if args.detect == 'Q':
        swirling = identification.q_criterion(a)
    elif args.detect == 'swirling':
        swirling = identification.calc_swirling(a)
    print(round(time.time() - lap,3), 'seconds')

    #---- PEAK DETECTION ----#
    print("Detecting peak of swirling strength")
    print("threshold=",args.threshold,"box size=",args.boxsize)

    peaks = detection.find_peaks(swirling, args.threshold, args.boxsize)

    print("Vortices found:",len(peaks[0]))

    #---- PEAKS DIRECTION OF ROTATION ----#
    dirL, dirR = detection.direction_rotation(vorticity,peaks)

    #---- SAVING OUTPUT FILE ----#
    if args.outfilename == None:
        pass
    else:
        print("saving file",args.outfilename)
    

  
    #---- PLOTTING ----#
    if args.plot_x == 'detect':
        plot.plot_detection(dirL,dirR,swirling)
    elif args.plot_x == 'fields':
        plot.plot_fields(a,vorticity)
    elif args.plot_x == 'quiver':
        for i in range(len(peaks[0])):
            xCenter = peaks[0][i]
            yCenter = peaks[1][i]
            dist = 10
            if (xCenter > dist) and (yCenter > dist):
                print('x1:',xCenter,'x2:',yCenter, 'swirl:',peaks[2][i])
                totalvel = np.sqrt(a.u**2+a.v**2)
                
                plot.plot_quiver(a, xCenter, yCenter, dist, swirling)
    else:
        print('no plot')


  
