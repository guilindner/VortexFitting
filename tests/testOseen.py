#!/usr/bin/env/ python3
"""
vortex detection tool, by Guilherme Lindner, 2017-04\n
This program load NetCDF files from DNS simulations  or PIV experiments
and detect the vortices and apply a fitting to them.
"""

import argparse
import numpy as np
import sys

sys.path.insert(1, '../vortexfitting')

import fitting  # noqa: E402
import schemes  # noqa: E402
import detection  # noqa: E402
from classes import VelocityField  # noqa: E402

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Optional app description',
                                     formatter_class = argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', dest = 'infilename',
                        default = '../data/example_dataHIT.nc',
                        help = 'input NetCDF file', metavar = 'FILE')
    args = parser.parse_args()

    print('Some tests for the Lamb-Oseen model')
    print(args.infilename)

    vfield = VelocityField(args.infilename, 0, '/', 'dns')


    def test_oseen(core_radius, gamma, dist, xdrift, ydrift, u_advection, v_advection):
        print('core_radius:', core_radius, 'Gamma', gamma, 'xdrift', xdrift,
              'ydrift', ydrift, 'u_advection', u_advection, 'v_advection', v_advection)
        model = [[], [], [], [], [], []]
        model[0] = core_radius
        model[1] = gamma
        core_radius_ori = model[0]
        gamma_ori = model[1]
        x_index = np.linspace(-1, 1, dist)
        y_index = np.linspace(-1, 1, dist)
        x_index, y_index = np.meshgrid(x_index, y_index)
        x_real = 0.0
        y_real = 0.0
        model[4] = u_advection
        model[5] = v_advection
        u_data, v_data = fitting.velocity_model(core_radius, gamma, x_real, y_real,
                                                u_advection, v_advection, x_index + xdrift, y_index + ydrift)
        u_data = u_data + u_advection
        v_data = v_data + v_advection
        # NOISE
        u_data = np.random.normal(u_data, 0.3)
        v_data = np.random.normal(v_data, 0.3)
        model = fitting.fit(core_radius, gamma, x_index, y_index, x_real, y_real, u_data, v_data, u_advection,
                            v_advection, 0)
        print('core_radius:', model[0], 'error(%):', (1 - (model[0]) / core_radius_ori) * 100)
        print('gamma:', model[1], 'error(%):', (1 - (model[1]) / gamma_ori) * 100)
        print('x_real:', model[2])
        print('y_real:', model[3])
        u_model, v_model = fitting.velocity_model(model[0], model[1], model[2], model[3],
                                                  model[4], model[5], x_index, y_index)
        corr = fitting.correlation_coef(u_data, v_data, u_model, v_model)
        print('correlation:', corr)
        print('---')
        fitting.plot_fit(x_index, y_index, u_data, v_data, u_model, v_model, model[2], model[3], model[0], model[1],
                         model[4], model[5], corr, 0, 0, '.', 0)


    test_oseen(0.2, 10, 10, 0.0, 0.0, 0.01, 0.01)
    test_oseen(0.2, 10, 10, 0.2, 0.2, 0.01, 0.01)
    test_oseen(0.9, 40, 10, 0.0, 0.0, 0.01, 0.01)
    test_oseen(0.9, 40, 10, 0.2, 0.2, 0.01, 0.01)
