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
import schemes
import detection
from classes import VelocityField

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Optional app description',
                                     formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-i', '--input', dest='infilename',
                        default='../data/test_data.nc',
                        help='input NetCDF file', metavar='FILE')
                        
    
    args = parser.parse_args()
    
    a = VelocityField(args.infilename)

    
    
    def test_oseen(coreR, gamma, dist):
        print('coreR:',coreR,'Gamma',gamma,'dist',dist)
        X = np.linspace(-1,1,dist)
        Y = np.linspace(-1,1,dist)
        a.dx = np.zeros(dist)
        a.dy = np.zeros(dist)
        X, Y = np.meshgrid(X,Y)
        xCenter = 0
        yCenter = 0
        xdrift = 0
        ydrift = 0
        coreRori = coreR
        gammaori = gamma
        Uw, Vw = fitting.velocity_model(a, X, Y,xCenter+xdrift,yCenter+ydrift, gamma, coreR)
        u_conv = 0.0 #flipped with v, fix later
        v_conv = 0.0
        Uw = Uw + u_conv
        Vw = Vw + v_conv
        coreR, gamma = fitting.fit(a, X, Y, xCenter, yCenter, Uw, Vw, u_conv, v_conv, gamma)
        print('coreR:',coreR,'error(%):',(1-(coreR)/coreRori)*100)
        print('gamma:',gamma,'error(%):',(1-(gamma)/gammaori)*100)
        #print('xCenter:', fxCenter)
        #print('yCenter:',fyCenter)
        uMod, vMod = fitting.velocity_model(a, X, Y,xCenter, yCenter, gamma, coreR)
        corr = fitting.correlation_coef(Uw,Vw,uMod,vMod)
        print('correlation:',corr)
        print('---')
        #plot.plot_corr(X, Y, Uw, Vw, uMod, vMod, coreR, corr)
  
    test_oseen(0.1,20,20)
    test_oseen(0.2,20,20)
    test_oseen(0.3,20,20)
    test_oseen(0.4,20,20)
    test_oseen(0.2,10,20)
    test_oseen(0.2,20,20)
    test_oseen(0.2,30,20)
    test_oseen(0.2,40,20)
