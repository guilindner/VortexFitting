import numpy as np
from scipy.signal import correlate2d
from scipy import optimize
from scipy.stats import pearsonr

import tools
import plot

def correlation_coef(Uw,Vw,u,v):
    """Calculates the correlation coefficient between two 2D arrays
    
    :param Uw: velocity u from the data at the proposed window
    :param Vw: velocity v from the data at the proposed window
    :param u: velocity u from the calculated model
    :param v: velocity v from the calculated model
    :type Uw: float
    :type Vw: float
    :type u: float
    :type v: float
    :returns: corr
    :rtype: float    
    """
    Uw = Uw.ravel()
    Vw = Vw.ravel()
    u = u.ravel()
    v = v.ravel()
    N = Uw.size
    prod_PIV_mod = 0.0
    prod_PIV = 0.0
    prod_mod = 0.0
    for i in range(N):
        prod_PIV_mod += (Uw[i]*u[i]+Vw[i]*v[i])/N
        prod_PIV += (u[i]*u[i]+v[i]*v[i])/N
        prod_mod += (Uw[i]*Uw[i]+Vw[i]*Vw[i])/N
    upper = prod_PIV_mod
    lower = np.sqrt(prod_PIV)*np.sqrt(prod_mod)
    corr = np.sqrt(upper/lower)

    return corr

def velocity_model(coreR, gamma, fxCenter,fyCenter, u_conv, v_conv,x,y):
    """Generates the Lamb-Oseen vortex velocity array
    
    :param coreR: core radius of the vortex
    :param gamma: circulation contained in the vortex
    :param fxCenter: relative x position of the vortex center
    :param fyCenter: relative y position of the vortex center
    :param u_conv: u convective velocity at the center
    :param v_conv: v convective velocity at the center
    :type coreR: float
    :type gamma: float
    :type fxCenter: float
    :type fyCenter: float
    :type u_conv: float
    :type v_conv: float
    :returns: velx, vely
    :rtype: float
    """
    r = np.hypot(x-fxCenter, y-fyCenter)
    vel = (gamma/(2 * np.pi * r)) * (1 - np.exp(-(r**2)/(coreR)**2))
    vel = np.nan_to_num(vel)
    velx = u_conv -vel*(y-fyCenter)/r
    vely = v_conv +vel*(x-fxCenter)/r
    velx = np.nan_to_num(velx)
    vely = np.nan_to_num(vely)
    return velx, vely

def get_vortices(a,peaks,vorticity):
    """General routine to check if the detected vortex is a real vortex
    
    :param a: data from the input file
    :param peaks: list of vortices
    :param vorticity: calculated field
    :type a: int, float,...
    :type peaks: list
    :type vorticity: array
    :returns: vortices
    :rtype: list
    """
    b = list()
    vortices = list()
    j = 0
    for i in range(len(peaks[0])):
        xCenter = peaks[1][i]
        yCenter = peaks[0][i]
        print(i," Processing detected swirling at (x,y)",xCenter,yCenter)
        coreR = 2*(a.dx[5]-a.dx[4]) #ugly change someday
        gamma = vorticity[yCenter,xCenter]*np.pi*coreR**2
        b = full_fit(coreR, gamma, a, xCenter, yCenter)
        X, Y, Uw, Vw = tools.window(a,b[2],b[3],b[6])
        uMod, vMod = velocity_model(b[0], b[1], b[2], b[3], b[4], b[5],X,Y)
        corr = correlation_coef(Uw,Vw,uMod,vMod)
        #print(corr)
        if (corr > 0.90):
            print("Accepted! corr = %s (vortex %s)" %(corr,j))
            vortices.append([b[0],b[1],b[2],b[3],b[4],b[5],b[6],corr])
            j += 1
    return vortices

def full_fit(coreR, gamma, a, xCenter, yCenter):
    """Full fitting procedure
    
    :param coreR: core radius of the vortex
    :param gamma: circulation contained in the vortex
    :param a: data from the input file
    :param xCenter: x index of the vortex center
    :param yCenter: y index of the vortex center
    :type coreR: float
    :type gamma: float
    :type a: class
    :type xCenter: int
    :type yCenter: int
    :returns: fitted[i],dist
    :rtype: list
    """

    fitted = [[],[],[],[],[],[]]
    fitted[0] = coreR
    fitted[1] = gamma
    fitted[2] = a.dx[xCenter]
    fitted[3] = a.dy[yCenter]
    dx = a.dx[5]-a.dx[4] #ugly
    dy = a.dy[5]-a.dy[4]
    corr = 0.0

    for i in range(10):
        xCenter = int(round(fitted[2]/dx))
        yCenter = int(round(fitted[3]/dy))
        if xCenter >= a.u.shape[1]:
            xCenter = a.u.shape[1]-1
        if yCenter >= a.v.shape[0]:
            yCenter = a.v.shape[0]-1
        r1 = fitted[0]
        x1 = fitted[2]
        y1 = fitted[3]
        dist = int(round(fitted[0]/dx,0)) + 1
        fitted[4] = a.u[yCenter, xCenter] #u_conv
        fitted[5] = a.v[yCenter, xCenter] #v_conv
        X, Y, Uw, Vw = tools.window(a,xCenter,yCenter,dist)
        fitted = fit(fitted[0], fitted[1], X, Y, fitted[2], fitted[3], Uw, Vw, fitted[4], fitted[5],i)
        if i > 0:
            if abs(fitted[0]/r1 -1) < 0.1:
                if (abs((fitted[2]/x1 -1)) < 0.1) or (abs((fitted[3]/y1 -1)) < 0.01):
                    #print("break")
                    break
            if (abs((fitted[2]-x1)) > dist*dx) or (abs((fitted[3]-y1)) > dist*dy):
                corr = 0.0
                break 
    return fitted[0], fitted[1], fitted[2], fitted[3], fitted[4], fitted[5], dist
def fit(coreR, gamma, x, y, fxCenter, fyCenter, Uw, Vw, u_conv, v_conv,i):
    """
    Fitting  of the Lamb-Oseen Vortex
    
    :param coreR: core radius of the vortex
    :param gamma: circulation contained in the vortex
    :param x: x position
    :param y: y position
    :param fxCenter: x position of the vortex center
    :param fyCenter: y position of the vortex center
    :type coreR: float
    :type gamma: float
    :type x: float
    :type y: float
    :type fxCenter: float
    :type fyCenter: float
    :returns: sol.x
    :rtype: float 
    """

    x = x.ravel()
    y = y.ravel()
    Uw = Uw.ravel()
    Vw = Vw.ravel()
    dx = x[1]-x[0]
    dy = dx
    def fun(fitted): #fitted[0]=coreR, fitted[1]=gamma, fitted[2]=fxCenter, fitted[3]=fyCenter, fitted[4]=u_conv, fitted[5]=v_conv
        """
        Lamb-Ossen velocity model used for the nonlinear fitting
        """
        r = np.hypot(x-fitted[2], y-fitted[3])
        expr2 = np.exp(-r**2/fitted[0]**2)
        z = fitted[1]/(2*np.pi*r) * (1 - expr2)
        z = np.nan_to_num(z)
        zx = fitted[4] -z*(y-fitted[3])/r -Uw
        zy = fitted[5] +z*(x-fitted[2])/r -Vw
        zx = np.nan_to_num(zx)
        zy = np.nan_to_num(zy)
        zt = np.append(zx,zy)
        return zt
    if i > 0:
        m = 1.0
    else:
        m = 4.0
    bnds=([coreR-coreR*m,gamma-abs(gamma)*m/2,fxCenter-m*dx,fyCenter-m*dy,u_conv-abs(u_conv),v_conv-abs(v_conv)],
          [coreR+coreR*m,gamma+abs(gamma)*m/2,fxCenter+m*dx,fyCenter+m*dy,u_conv+abs(u_conv),v_conv+abs(v_conv)])
    sol = optimize.least_squares(fun, [coreR,gamma,fxCenter,fyCenter,u_conv,v_conv],bounds=bnds)     
    return sol.x
