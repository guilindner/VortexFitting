#!/usr/bin/env/ python
"""vortex detection tool, by Guilherme Lindner, 2017-04
This program load NetCDF files from DNS simulations  or PIV experiments
and detect the vortices and apply a fitting to them.
"""
import numpy as np

def calc_swirling(a):
    A = np.zeros((a.sizex*a.sizey,3,3))
    A = np.array([[a.derivative['dudx'].ravel(),a.derivative['dudy'].ravel(),
                a.derivative['dudz'].ravel()],[a.derivative['dvdx'].ravel(),
                a.derivative['dvdy'].ravel(),a.derivative['dvdz'].ravel()],
                [a.derivative['dwdx'].ravel(),a.derivative['dwdy'].ravel(),
                -a.derivative['dudx'].ravel()-a.derivative['dvdy'].ravel()]])
    A = A.transpose(2,1,0)
    eigenvalues = np.linalg.eigvals(A)
    swirling = np.max(eigenvalues.imag,axis=1).reshape(a.sizex,a.sizey)
    return swirling
