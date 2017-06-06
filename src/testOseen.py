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
#import detection
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
                        default=1.5, type=float,
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
    
    a = VelocityField(args.infilename,args.timestep)

    xCenter = 0
    yCenter = 0
    gamma = -100
    coreR = 0.9
    dist = 20

    #X, Y, Uw, Vw = tools.window(a,xCenter,yCenter,dist)
    X = np.linspace(-1,1,dist)
    Y = np.linspace(-1,1,dist)
    a.dx = np.zeros(dist)
    a.dy = np.zeros(dist)
    X, Y = np.meshgrid(X,Y)   
    uMod, vMod = fitting.velocity_model(a, X, Y,xCenter,yCenter, gamma, coreR)
    Uw, Vw = fitting.velocity_model(a, X, Y,xCenter,yCenter, gamma, coreR)
    u_conv = 0.0 #flipped with v, fix later
    v_conv = 0.0
    Uw = Uw + u_conv
    Vw = Vw + v_conv
    
    corr = fitting.correlation_coef(Uw,Vw,uMod,vMod)
    print('correlation before fit:',corr)
    coreR, gamma = fitting.fit(a, X, Y, xCenter, yCenter, Uw, Vw, u_conv, v_conv, gamma)
    print('coreR:',coreR)
    uMod, vMod = fitting.velocity_model(a, X, Y,xCenter,yCenter, gamma, coreR)
    corr = fitting.correlation_coef(Uw,Vw,uMod,vMod)
    print('correlation AFTER fit:',corr,'and gamma!',gamma)
    plot.plot_corr(X, Y, Uw, Vw, uMod, vMod, coreR, corr)
  
