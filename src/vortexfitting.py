#!/usr/bin/env/ python3
"""vortex detection tool, by Guilherme Lindner, 2017-04\n
This program load NetCDF files from DNS simulations  or PIV experiments
and detect the vortices and apply a fitting to them.
"""
import sys
import argparse
import time
import numpy as np

from classes import VelocityField
import tools
import fitting
import plot
import schemes
import detection

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Optional app description',
                                     formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-i', '--input', dest='infilename',
                        default='../data/test_data.nc',
                        help='input NetCDF file', metavar='FILE')
    
    parser.add_argument('-s', '--scheme', dest='scheme', type=int, default=2,
                        help='Scheme for differencing\n'
                             '2 = second order\n'
                             '22 = least-square filter'
                             '4 = fourth order')
    
    parser.add_argument('-T', '--time', dest='timestep', type=int,
                        default=0,
                        help='Timestep/Sample/Z position desired')
                        
    parser.add_argument('-d', '--detect', dest='detect',
                        default='swirling',
                        help='Detection method:\n'
                             'Q = Q criterion\n'
                             'delta = delta criterion\n'
                             'swirling = 2D Swirling Strength')
    
    parser.add_argument('-t', '--threshold', dest='threshold',
                        default=0.1, type=float,
                        help='Threshold for detection, integer')

    parser.add_argument('-b', '--boxsize', dest='boxsize',
                        default=6, type=int,
                        help='Box size for the detection')
    
    parser.add_argument('-f', '--flip', dest='flip',
                        default=False, type=bool,
                        help='Flip X and Y axis for plotting, 0 = False, 1 = True')
    
    parser.add_argument('-p', '--plot', dest='plot_x',
                        default='fit',
                        help='Plot on screen:\n'
                             'fit    = Detection and fitting, saves images (default)\n'
                             'detect = Possible vortices (no fitting)\n'
                             'fields = Velocity fields and vorticity\n')
    parser.add_argument('-xy', '--xy', nargs=2, dest='xy', default=[0,0])
    
    args = parser.parse_args()
    
    start = time.time()
    #---- LOAD DATA ----#
    print("Opening file:",args.infilename)

    #print("Sample target: (todo)", args.timestep)
    
    a = VelocityField(args.infilename,args.timestep)
    print("Samples:", a.samples)

    #---- DIFFERENCE APPROXIMATION ----# 
    lap = time.time()
    if args.scheme == 4:
        a.derivative = schemes.fourth_order_diff(a)
    elif args.scheme == 2:
        a.derivative = schemes.second_order_diff(a)
    elif args.scheme == 22:
        a.derivative = schemes.least_square_diff(a)
    else:
        print('No scheme', args.scheme, 'found. Exitting!')
        sys.exit()
    #print(round(time.time() - lap,3), 'seconds') 
    
    #---- VORTICITY ----#

    vorticity = a.derivative['dvdx'] - a.derivative['dudy']

    #---- METHOD FOR DETECTION OF VORTICES ----#
    lap = time.time()
    if args.detect == 'Q':
        swirling = detection.q_criterion(a)
    elif args.detect == 'swirling':
        swirling = detection.calc_swirling(a)
    elif args.detect == 'delta':
        swirling = detection.delta_criterion(a)
    #print(round(time.time() - lap,3), 'seconds')

    if a.norm == True:
        swirling = tools.normalize(swirling,a.normdir) #normalization

    #---- PEAK DETECTION ----#
    print("threshold=",args.threshold,"box size=",args.boxsize)

    peaks = tools.find_peaks(swirling, args.threshold, args.boxsize)

    print("Vortices found:",len(peaks[0]))
    #---- PEAKS DIRECTION OF ROTATION ----#
    dirL, dirR = tools.direction_rotation(vorticity,peaks)

    #---- MODEL FITTING ----# SEE IN PLOT
    vortices = list()
    if (args.plot_x != 'fit') or (args.xy != [0,0]):
        print("No fitting")
    else:
        vortices = fitting.get_vortices(a,peaks,vorticity)
        print('---- Accepted vortices ----')
        print(len(vortices))
    #print('xCenter, yCenter, gamma, core Radius, correlation, mesh distance')
    #for vortex in vortices:
    #    print(vortex)    
  
    #---- PLOTTING OPTIONS ----#
    if args.xy != [0,0]:
        x = int(args.xy[0])
        y = int(args.xy[1])
        swirlingw = swirling[y-10:y+10,x-10:x+10]
        X, Y, Uw, Vw = tools.window(a,x,y,10)
        plot.plot_quiver(X, Y, Uw, Vw, swirlingw)
    elif args.plot_x == 'detect':
        plot.plot_detect(dirL,dirR,swirling,args.flip)
    elif args.plot_x == 'fields':
        plot.plot_fields(a,vorticity)
    elif args.plot_x == 'fit':
        plot.plot_accepted(vortices,swirling)
        outfile = open('../results/vortices.dat','w')
        outfile.write('X Y gamma radius corr mesh x y u_c v_c \n')
        for line in vortices:
            outfile.write("%s %s %s %s %s %s %s %s %s %s \n" % line)
        for i in range(len(vortices)):
            print('x:',vortices[i][0],'y:',vortices[i][1], 'gamma:',vortices[i][2],
             'mesh',vortices[i][5], 'corr',vortices[i][4], 'coreR',vortices[i][3])
            X, Y, Uw, Vw = tools.window(a,vortices[i][0],vortices[i][1],vortices[i][5]*2)
            uMod, vMod = fitting.velocity_model(vortices[i][3], vortices[i][2],
             vortices[i][6], vortices[i][7], vortices[i][8], vortices[i][9], X, Y)
            corr = fitting.correlation_coef(Uw,Vw,uMod,vMod)
            plot.plot_fit(X, Y, Uw, Vw, uMod, vMod, vortices[i][6],
                       vortices[i][7], vortices[i][3], vortices[i][2], vortices[i][4],i)
        
    else:
        print('no plot')
