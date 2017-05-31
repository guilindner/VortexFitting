#!/usr/bin/env/ python3
"""vortex detection tool, by Guilherme Lindner, 2017-04\n
This program load NetCDF files from DNS simulations  or PIV experiments
and detect the vortices and apply a fitting to them.
"""
import sys
import argparse
import time
import numpy as np
import scipy
from scipy import optimize

import tools
import fitting
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
                        default=11., type=float,
                        help='Threshold for detection, integer')

    parser.add_argument('-b', '--boxsize', dest='boxsize',
                        default=6, type=int,
                        help='Box size for the detection')
    
    parser.add_argument('-p', '--plot', dest='plot_x',
                        default='',
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

    #---- MODEL FITTING ----# SEE IN PLOT

    dist = 5
    coreR = 0.3


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
            X, Y, Uw, Vw = tools.window(a,xCenter,yCenter,dist)
            swirlingw = swirling[xCenter-dist:xCenter+dist,yCenter-dist:yCenter+dist] #reuse window function?
            if (xCenter > dist) and (yCenter > dist):
                print('x1:',xCenter,'x2:',yCenter, 'swirl:',peaks[2][i])
                plot.plot_quiver(X, Y, Uw, Vw, swirlingw)
    elif args.plot_x == 'corr':
        for i in range(len(peaks[0])):
            xCenter = peaks[0][i]
            yCenter = peaks[1][i]
            gamma = vorticity[xCenter,yCenter]
            X, Y, Uw, Vw = tools.window(a,xCenter,yCenter,dist)
            uMod, vMod = fitting.velocity_model(a, X, Y,xCenter,yCenter, vorticity[xCenter,yCenter], coreR)
            if (xCenter > dist) and (yCenter > dist):
                print('x1:',xCenter,'x2:',yCenter, 'swirl:',peaks[2][i])
                corr = fitting.correlation_coef(Uw,Vw,uMod,vMod)
                print('Correlation R =',corr)
                if (corr > 0.75):
                    print('R > 0.75, it\'s a vortex')
                    print('uMod',uMod)
                    print('Uw',Uw)
                else:
                    print('not a vortex')
                plot.plot_corr(X, Y, Uw, Vw, uMod, vMod)
    elif args.plot_x == 'fit':
        for i in range(1):
            xCenter = peaks[0][i]
            yCenter = peaks[1][i]
            
            gamma = vorticity[xCenter,yCenter]
            X, Y, Uw, Vw = tools.window(a,xCenter,yCenter,dist)
            uMod, vMod = fitting.velocity_model(a, X, Y,xCenter,yCenter, vorticity[xCenter,yCenter], coreR)
            if (xCenter > dist) and (yCenter > dist):
                print('xCenter:',xCenter,'yCenter:',yCenter, 'vorticity:',gamma)
                corr = fitting.correlation_coef(Uw,Vw,uMod,vMod)
                print('Correlation R =',corr)
                
                if (corr > 0.75):
                    print('R > 0.75, it\'s a vortex')
                else:
                    print('not a vortex, initializing fitting')
                    fitx = fitting.super_fitx(X, Y, Uw, Vw, gamma)
                    fity = fitting.super_fity(X, Y, Uw, Vw, gamma)
                    #fit = optimize.root(fitting.funx, 2.0, jac=fitting.jacx, method='lm')
                    print(fitx,fity)
                    uMod, vMod = fitting.velocity_model(a, X, Y,xCenter,yCenter, vorticity[xCenter,yCenter], fity)
                    corr = fitting.correlation_coef(Uw,Vw,uMod,vMod)
                    #print(Uw)
                    #print(uMod)
                    print('New Correlation R =',corr)
                    plot.plot_corr(X, Y, Uw, Vw, uMod, vMod)
                    plt.show()
                #fit = optimize.root(fitting.model_oseen_x, 1.0, method='lm')
                #root = root(model_oseen(a, x, y,xCenter,yCenter, gamma, coreR),method='lm')
                #print(fit)
                #plot.plot_corr(X, Y, Uw, Vw, uMod, vMod)
                
    else:
        print('no plot')


  
