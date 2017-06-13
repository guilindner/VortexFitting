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

    
    
    def test_oseen(coreR, gamma, dist,xdrift,ydrift):
        print('|*|coreR:',coreR,'Gamma',gamma,'xdrift',xdrift,'ydrift',ydrift,'|*|')
        X = np.linspace(-1,1,dist)
        Y = np.linspace(-1,1,dist)
        a.dx = np.zeros(dist)
        a.dy = np.zeros(dist)
        X, Y = np.meshgrid(X,Y)
        fxCenter = 0.0
        fyCenter = 0.0
        coreRori = coreR
        gammaori = gamma
        u_conv = 0.2 #flipped with v, fix later
        v_conv = 0.0
        Uw, Vw = fitting.velocity_model(u_conv, v_conv, X+xdrift, Y+ydrift,fxCenter,fyCenter, gamma, coreR)
        Uw = Uw + u_conv
        Vw = Vw + v_conv
        # NOISE
        Uw = np.random.normal(Uw,np.max(Uw)*0.1)
        Vw = np.random.normal(Vw,np.max(Vw)*0.1)
        coreR, gamma, fxCenter, fyCenter, u_conv, v_conv = fitting.fit(a, X, Y, fxCenter, fyCenter, Uw, Vw, u_conv, v_conv, gamma)
        print('coreR:',coreR,'error(%):',(1-(coreR)/coreRori)*100)
        print('gamma:',gamma,'error(%):',(1-(gamma)/gammaori)*100)
        print('xCenter:',fxCenter)
        print('yCenter:',fyCenter)
        print('u_conv:',u_conv)
        print('v_conv:',v_conv)
        #print('xCenter:', fxCenter)
        #print('yCenter:',fyCenter)
        uMod, vMod = fitting.velocity_model(u_conv, v_conv, X, Y, fxCenter, fyCenter, gamma, coreR)
        corr = fitting.correlation_coef(Uw,Vw,uMod,vMod)
        print('correlation:',corr)
        print('---')
        plot.plot_corr(X, Y, Uw, Vw, uMod, vMod, coreR, corr)
  
    test_oseen(0.2,10,10,0.0,0.0)
    test_oseen(0.2,-10,10,0.0,0.0)
    test_oseen(0.2,1,10,0.2,0.2)
    test_oseen(0.9,40,10,0.0,0.0)
    test_oseen(0.9,40,10,0.2,0.2)
