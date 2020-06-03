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

    vfield = VelocityField(args.infilename,0,'/','dns')

    def test_oseen(coreR, gamma, dist,xdrift,ydrift,u_conv,v_conv):
        print('coreR:',coreR,'Gamma',gamma,'xdrift',xdrift,
              'ydrift',ydrift,'u_conv',u_conv,'v_conv',v_conv)
        model = [[],[],[],[],[],[]]
        model[0] = coreR
        model[1] = gamma
        coreRori = model[0]
        gammaori = model[1]
        x_index = np.linspace(-1,1,dist)
        y_index = np.linspace(-1,1,dist)
        x_index, y_index = np.meshgrid(x_index, y_index)
        x_real = 0.0
        y_real = 0.0
        model[4] = u_conv
        model[5] = v_conv
        u_data, v_data = fitting.velocity_model(coreR, gamma, x_real, y_real,
                                                u_conv, v_conv, x_index+xdrift, y_index+ydrift)
        u_data = u_data + u_conv
        v_data = v_data + v_conv
        # NOISE
        u_data = np.random.normal(u_data,0.3)
        v_data = np.random.normal(v_data,0.3)
        model = fitting.fit(coreR, gamma, x_index, y_index, x_real, y_real, u_data, v_data, u_conv, v_conv,0)
        print('coreR:',model[0],'error(%):',(1-(model[0])/coreRori)*100)
        print('gamma:',model[1],'error(%):',(1-(model[1])/gammaori)*100)
        print('x_real:',model[2])
        print('y_real:',model[3])
        u_model, v_model = fitting.velocity_model(model[0], model[1], model[2], model[3],
                                                  model[4],model[5], x_index, y_index)
        corr = fitting.correlation_coef(u_data,v_data,u_model,v_model)
        print('correlation:',corr)
        print('---')
        plot.plot_fit(x_index, y_index, u_data, v_data, u_model, v_model, model[2], model[3], model[0], model[1], model[4],model[5], corr,0,0,'.',0)

    test_oseen(0.2,10,10,0.0,0.0,0.01,0.01)
    test_oseen(0.2,10,10,0.2,0.2,0.01,0.01)
    test_oseen(0.9,40,10,0.0,0.0,0.01,0.01)
    test_oseen(0.9,40,10,0.2,0.2,0.01,0.01)
