import numpy as np

def second_order_diff(a):
    """Second order accurate finite difference scheme.

    .. note:: Scheme 1:0:-1

    :params a: 2D array of the velocity field, containing u and v
    :type a: float
    :returns: derivative
    :rtype: float
    """
    print("Difference scheme: Second Order Scheme")
    x, y = np.meshgrid(a.dx,a.dy)
    dx = a.dx[1]-a.dx[0] #only for homogeneous mesh
    dy = a.dy[1]-a.dy[0] #only for homogeneous mesh
    a.derivative['dudy'], a.derivative['dudx'] = np.gradient(a.u,dx)
    a.derivative['dvdy'], a.derivative['dvdx'] = np.gradient(a.v,dy)
    return a.derivative

def least_square_diff(a): #there is a problem on the boundary
    """Fourth order accurate finite difference scheme.

    Least-square filter (Raffel 1998)

    .. note:: Scheme -2:-1:0:1:2

    :params a: 2D array of the velocity field, containing u and v
    :type a: float
    :returns: derivative
    :rtype: float
    """
    print("Difference scheme: least-square filter")
    dx = a.dx[1]-a.dx[0] #only for homogeneous mesh
    dy = a.dy[1]-a.dy[0] #only for homogeneous mesh
    ### INVERT AXIS!!! ###
    a.derivative['dudx'][:,2:-2] = (-2*a.u[:, 0:-4] - a.u[:,1:-3]+ a.u[:, 3:-1] + 2*a.u[:,4:])/(10*dy)
    a.derivative['dudy'][2:-2,:] = (-2*a.u[0:-4,:] - a.u[1:-3,:] + a.u[3:-1,:] + 2*a.u[4:,:])/(10*dx)
    a.derivative['dvdx'][:,2:-2] = (-2*a.v[:, 0:-4] - a.v[:,1:-3]+ a.v[:, 3:-1] + 2*a.v[:,4:])/(10*dy)
    a.derivative['dvdy'][2:-2,:] = (-2*a.v[0:-4,:] - a.v[1:-3,:] + a.v[3:-1,:] + 2*a.v[4:,:])/(10*dx)

    return a.derivative


def fourth_order_diff(a):
    """Fourth order accurate finite difference scheme.

    .. note:: Scheme: 1:-8:0:8:-1

    :params a: 2D array of the velocity field, containing u and v
    :type a: float
    :returns: derivative
    :rtype: float
    """
    print("Difference scheme: Fourth Order Scheme")
    dx = a.dx[1]-a.dx[0] #only for homogeneous mesh
    dy = a.dy[1]-a.dy[0] #only for homogeneous mesh
    a.derivative['dudx'][:,2:-2] = (a.u[:, 0:-4] -8*a.u[:,1:-3]+ 8*a.u[:, 3:-1] - a.u[:,4:])/(12*dy)
    a.derivative['dudy'][2:-2,:] = (a.u[0:-4,:] -8*a.u[1:-3,:] + 8*a.u[3:-1,:] - a.u[4:,:])/(12*dx)
    a.derivative['dvdx'][:,2:-2] = (a.v[:, 0:-4] -8*a.v[:,1:-3]+ 8*a.v[:, 3:-1] - a.v[:,4:])/(12*dy)
    a.derivative['dvdy'][2:-2,:] = (a.v[0:-4,:] -8*a.v[1:-3,:] + 8*a.v[3:-1,:] - a.v[4:,:])/(12*dx)

    return a.derivative
