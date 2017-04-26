"""Finite differences schemes
"""

import numpy as np

def secondOrderDiff(a):
    """Second order accurate finite difference scheme.
    
    Scheme 1:0:-1
    
    Args:
        a: 2D array of the velocity field, containing u and v
    
    Returns:
        All the spatial derivatives, dudx, dudy, dvdx and dvdy
        For example:
        a.derivative['dudx']
    """
    print("Beginning differenciation with Second Order Scheme")
    a.derivative['dudx'][0,0] = (a.u[1,0] - a.u[0,0])/(a.dx[1]-a.dx[0])
    a.derivative['dudx'][-1,0] = (a.u[-1,0] - a.u[-2,0])/(a.dx[-1]-a.dx[-2])
#    a.derivative['dudy'][0,0] = (a.u[1,0] - a.u[0,0])/(a.dx[1]-a.dx[0])
    a.derivative['dudx'][1:-1,1:-1] = (a.u[2:, 1:-1] - a.u[:-2,1:-1])/(2*(a.dx[1:-1,np.newaxis]-a.dx[0:-2,np.newaxis]))
    a.derivative['dudy'][1:-1,1:-1] = (a.u[1:-1, 2:] - a.u[1:-1,:-2])/(2*(a.dy[1:-1]-a.dy[0:-2]))
    a.derivative['dvdx'][1:-1,1:-1] = (a.v[2:, 1:-1] - a.v[:-2,1:-1])/(2*(a.dx[1:-1,np.newaxis]-a.dx[0:-2,np.newaxis]))
    a.derivative['dvdy'][1:-1,1:-1] = (a.v[1:-1, 2:] - a.v[1:-1,:-2])/(2*(a.dy[1:-1]-a.dy[0:-2]))
    return a.derivative
    
    
def SLOWsecondOrderDiff(a): # 1200 times slower
    """Second order accurate finite difference scheme"""
    print("Beginning differenciation with Second Order Scheme")
    x = a.u.shape[0]-1
    y = a.u.shape[1]-1

    for i in range(1,x-1):
        for j in range(1,y-1):
            a.derivative['dudx'][i,j] = (a.u[i+1, j] - a.u[i-1,j])/(2*(a.dx[i]-a.dx[i-1]))
            a.derivative['dudy'][i,j] = (a.u[i, j+1] - a.u[i,j-1])/(2*(a.dy[j]-a.dy[j-1]))
            a.derivative['dvdx'][i,j] = (a.v[i+1, j] - a.v[i-1,j])/(2*(a.dx[i]-a.dx[i-1]))
            a.derivative['dvdy'][i,j] = (a.v[i, j+1] - a.v[i,j-1])/(2*(a.dy[j]-a.dy[j-1]))
    return a.derivative

def fourthOrderDiff(a):
    """Fourth order accurate finite difference scheme.
    
    Scheme: -1:8:0:-8:1
    
    Args:
        a: 2D array of the velocity field, containing u and v
    
    Returns:
        All the spatial derivatives, dudx, dudy, dvdx and dvdy
        For example:
        a.derivative['dudx']
    """
    print("Beginning differenciation with Fourth Order Scheme")
    #dont work for non-uniform mesh
    x = a.u.shape[0]-1
    y = a.u.shape[1]-1
    a.derivative['dudx'][1,1] = (a.u[4,1] -6*a.u[3,1] + 18*a.u[2,1] -10*a.u[1,1] -3*a.u[0,1])/(12*(a.dx[1]-a.dx[0]))
    a.derivative['dudx'][x-1,1] = (3*a.u[x,1] +10*a.u[x-1,1] -18*a.u[x-2,1] +6*a.u[x-3,1] -1*a.u[x-4,1])/(12*(a.dx[1]-a.dx[0]))
    a.derivative['dvdx'][1,1] = (a.v[4,1] -6*a.v[3,1] + 18*a.v[2,1] -10*a.v[1,1] -3*a.v[0,1])/(12*(a.dx[1]-a.dx[0]))
    a.derivative['dvdx'][x-1,1] = (3*a.v[x,1] +10*a.v[x-1,1] -18*a.v[x-2,1] +6*a.v[x-3,1] -1*a.v[x-4,1])/(12*(a.dx[1]-a.dx[0]))
    a.derivative['dudy'][1,1] = (a.u[1,4] -6*a.u[1,3] + 18*a.u[1,2] -10*a.u[1,1] -3*a.u[1,0])/(12*(a.dy[1]-a.dy[0]))
    a.derivative['dudy'][1,y-1] = (3*a.u[1,y] +10*a.u[1,y-1] -18*a.u[1,y-2] +6*a.u[1,y-3] -1*a.u[1,y-4])/(12*(a.dy[1]-a.dy[0]))
    a.derivative['dvdy'][1,1] = (a.v[1,4] -6*a.v[1,3] + 18*a.v[1,2] -10*a.v[1,1] -3*a.v[1,0])/(12*(a.dx[1]-a.dx[0]))
    a.derivative['dvdy'][1,y-1] = (3*a.v[1,y] +10*a.v[1,y-1] -18*a.v[1,y-2] +6*a.v[1,y-3] -1*a.v[1,y-4])/(12*(a.dy[1]-a.dy[0]))
    for i in range(2,x-2):
        for j in range(2,y-2):
            a.derivative['dudx'][i,j] = (-a.u[i+2,j] + 8*a.u[i+1, j] - 8*a.u[i-1,j] + a.u[i-2,j])/(12*(a.dx[i]-a.dx[i-1]))
            a.derivative['dudy'][i,j] = (-a.u[i,j+2] + 8*a.u[i, j+1] - 8*a.u[i,j-1] + a.u[i,j-2])/(12*(a.dy[i]-a.dy[i-1]))
            a.derivative['dvdx'][i,j] = (-a.v[i+2,j] + 8*a.v[i+1, j] - 8*a.v[i-1,j] + a.v[i-2,j])/(12*(a.dx[i]-a.dx[i-1]))
            a.derivative['dvdy'][i,j] = (-a.v[i,j+2] + 8*a.v[i, j+1] - 8*a.v[i,j-1] + a.v[i,j-2])/(12*(a.dy[i]-a.dy[i-1]))
    return a.derivative
    
#    a.derivative['dudx'][2:-2,2:-2] = (-a.u[4:,2:-2] + 8*a.u[2:, 2:-2] - 8*a.u[1:-3,] + a.u[i-2,j])/(12*(a.dx[i]-a.dx[i-1]))

#    a.derivative['dudx'][i,j] = (a.u[i+1, j] - a.u[i-1,j])/(2*(a.dx[i]-a.dx[i-1]))
#    a.derivative['dudx'][1:-1,1:-1] = (a.u[2:, 1:-1] - a.u[:-2,1:-1])/(2*(a.dx[1:-1,np.newaxis]-a.dx[0:-2,np.newaxis]))
