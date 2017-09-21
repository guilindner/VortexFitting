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
                        default='../data/test_dataHIT.nc',
                        help='input NetCDF file', metavar='FILE')

    args = parser.parse_args()

    a = VelocityField(args.infilename)

    def test_oseen(coreR, gamma, dist,xdrift,ydrift,u_conv,v_conv):
        print('coreR:',coreR,'Gamma',gamma,'xdrift',xdrift,
              'ydrift',ydrift,'u_conv',u_conv,'v_conv',v_conv)
        model = [[],[],[],[],[],[]]
        model[0] = coreR
        model[1] = gamma
        coreRori = model[0]
        gammaori = model[1]
        X = np.linspace(-1,1,dist)
        Y = np.linspace(-1,1,dist)
        X, Y = np.meshgrid(X,Y)
        fxCenter = 0.0
        fyCenter = 0.0
        model[4] = u_conv
        model[5] = v_conv
        Uw, Vw = fitting.velocity_model(coreR, gamma, fxCenter, fyCenter, u_conv, v_conv, X+xdrift, Y+ydrift)
        Uw = Uw + u_conv
        Vw = Vw + v_conv
        # NOISE
        Uw = np.random.normal(Uw,0.3)
        Vw = np.random.normal(Vw,0.3)
        model = fitting.fit(coreR, gamma, X, Y, fxCenter, fyCenter, Uw, Vw, u_conv, v_conv,0)
        print('coreR:',model[0],'error(%):',(1-(model[0])/coreRori)*100)
        print('gamma:',model[1],'error(%):',(1-(model[1])/gammaori)*100)
        print('fxCenter:',model[2])
        print('fyCenter:',model[3])
        uMod, vMod = fitting.velocity_model(model[0], model[1], model[2], model[3],model[4],model[5],X,Y)
        corr = fitting.correlation_coef(Uw,Vw,uMod,vMod)
        print('correlation:',corr)
        print('---')
        plot.plot_fit_test(X, Y, Uw, Vw, uMod, vMod, model[2], model[3], model[0], model[1], model[4],model[5], corr)

    test_oseen(0.2,10,10,0.0,0.0,0.01,0.01)
    test_oseen(0.2,10,10,0.2,0.2,0.01,0.01)
    test_oseen(0.9,40,10,0.0,0.0,0.01,0.01)
    test_oseen(0.9,40,10,0.2,0.2,0.01,0.01)
