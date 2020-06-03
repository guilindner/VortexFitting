import numpy as np

def second_order_diff(vfield):
    """Second order accurate finite difference scheme.

    .. note:: Scheme 1:0:-1

    :params vfield: 2D array of the velocity field, containing u and v
    :type vfield: float
    :returns: derivative
    :rtype: float
    """
    print("Difference scheme: Second Order Scheme")
    x, y = np.meshgrid(vfield.dx,vfield.dy)
    dx = vfield.dx[1]-vfield.dx[0] #only for homogeneous mesh
    dy = vfield.dy[1]-vfield.dy[0] #only for homogeneous mesh
    vfield.derivative['dudy'], vfield.derivative['dudx'] = np.gradient(vfield.u,dx)
    vfield.derivative['dvdy'], vfield.derivative['dvdx'] = np.gradient(vfield.v,dy)
    return vfield.derivative

def least_square_diff(vfield): #may have a problem on the boundary
    """Fourth order accurate finite difference scheme.

    Least-square filter (Raffel 1998)

    .. note:: Scheme -2:-1:0:1:2

    :params vfield: 2D array of the velocity field, containing u and v
    :type vfield: float
    :returns: derivative
    :rtype: float
    """
    print("Difference scheme: least-square filter")
    dx = vfield.dx[1]-vfield.dx[0] #only for homogeneous mesh
    dy = vfield.dy[1]-vfield.dy[0] #only for homogeneous mesh
    ### INVERT AXIS!!! ###
    vfield.derivative['dudx'][:,2:-2] = (-2*vfield.u[:, 0:-4] - vfield.u[:,1:-3]+ vfield.u[:, 3:-1] + 2*vfield.u[:,4:])/(10*dy)
    vfield.derivative['dudy'][2:-2,:] = (-2*vfield.u[0:-4,:] - vfield.u[1:-3,:] + vfield.u[3:-1,:] + 2*vfield.u[4:,:])/(10*dx)
    vfield.derivative['dvdx'][:,2:-2] = (-2*vfield.v[:, 0:-4] - vfield.v[:,1:-3]+ vfield.v[:, 3:-1] + 2*vfield.v[:,4:])/(10*dy)
    vfield.derivative['dvdy'][2:-2,:] = (-2*vfield.v[0:-4,:] - vfield.v[1:-3,:] + vfield.v[3:-1,:] + 2*vfield.v[4:,:])/(10*dx)

    return vfield.derivative


def fourth_order_diff(vfield):
    """Fourth order accurate finite difference scheme.

    .. note:: Scheme: 1:-8:0:8:-1

    :params vfield: 2D array of the velocity field, containing u and v
    :type vfield: float
    :returns: derivative
    :rtype: float
    """
    print("Difference scheme: Fourth Order Scheme")
    dx = vfield.dx[1]-vfield.dx[0] #only for homogeneous mesh
    dy = vfield.dy[1]-vfield.dy[0] #only for homogeneous mesh
    vfield.derivative['dudx'][:,2:-2] = (vfield.u[:, 0:-4] -8*vfield.u[:,1:-3]+ 8*vfield.u[:, 3:-1] - vfield.u[:,4:])/(12*dy)
    vfield.derivative['dudy'][2:-2,:] = (vfield.u[0:-4,:] -8*vfield.u[1:-3,:] + 8*vfield.u[3:-1,:] - vfield.u[4:,:])/(12*dx)
    vfield.derivative['dvdx'][:,2:-2] = (vfield.v[:, 0:-4] -8*vfield.v[:,1:-3]+ 8*vfield.v[:, 3:-1] - vfield.v[:,4:])/(12*dy)
    vfield.derivative['dvdy'][2:-2,:] = (vfield.v[0:-4,:] -8*vfield.v[1:-3,:] + 8*vfield.v[3:-1,:] - vfield.v[4:,:])/(12*dx)

    return vfield.derivative
