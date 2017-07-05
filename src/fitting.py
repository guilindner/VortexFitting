import numpy as np
from scipy.signal import correlate2d
from scipy import optimize
from scipy.stats import pearsonr

import tools
import plot

def correlation_coef(Uw,Vw,u,v):
    Rx = pearsonr(Uw.ravel(),u.ravel())
    Ry = pearsonr(Vw.ravel(),v.ravel())
    R = Rx[0]*Ry[0]
    
    return R

def velocity_model(coreR, gamma, fxCenter,fyCenter, u_conv, v_conv,x,y):
    r = np.hypot(x-fxCenter, y-fyCenter)
    vel = (gamma/(2 * np.pi * r)) * (1 - np.exp(-(r**2)/(coreR)**2))
    vel = np.nan_to_num(vel)
    velx = (vel + u_conv)*(-x+fxCenter)
    vely = (vel + v_conv)*(y-fyCenter)
    return velx, vely

def temporary(a,peaks,vorticity):
    b = list()
    vortices = list()
    for i in range(len(peaks[0])):
        print("Processing Vortex:",i)
        xCenter = peaks[0][i]
        yCenter = peaks[1][i]
        if (len(a.dx)-10 > xCenter > 10) and (len(a.dy)-10 > yCenter > 10): #skip near wall
            gamma = vorticity[xCenter,yCenter]
            coreR = 4*(a.dx[xCenter+1]-a.dx[xCenter])
            b = full_fit(coreR, gamma, a, xCenter, yCenter)
            if (b[4] > 0.75):
                print("Accepted!")
                vortices.append(b)
    return vortices

def full_fit(coreR, gamma, a, xCenter, yCenter):
    model = [[],[],[],[],[],[]]
    model[0] = coreR
    model[1] = gamma
    model[2] = a.dx[xCenter]
    model[3] = a.dy[yCenter]
    dx = a.dx[xCenter+1]-a.dx[xCenter]
    dy = a.dy[yCenter+1]-a.dy[yCenter]
    dist = int(round(model[0]/dx,0)) + 1
    u_conv = a.u[xCenter, yCenter]
    v_conv = a.v[xCenter, yCenter]
    X, Y, Uw, Vw = tools.window(a,xCenter,yCenter,dist)
    model = fit(model[0], model[1], X, Y, model[2], model[3], Uw, Vw, u_conv, v_conv)
    uMod, vMod = velocity_model(model[0], model[1], model[2], model[3], u_conv, v_conv,X,Y)
    corr = correlation_coef(Uw,Vw,uMod,vMod)
    #plot.plot_corr(X, Y, Uw, Vw, uMod, vMod, model[0], corr)
    
    if (xCenter > len(a.u[0])):
        xCenter = len(a.u[0])
    if (yCenter > len(a.v[0])):
        yCenter = len(a.v[0])        

    if (corr > 0.75):
        plot.plot_corr(X, Y, Uw, Vw, uMod, vMod, model[0], corr)
        dist = int(round(2*model[0]/dx,0))
        u_conv = a.u[xCenter, yCenter]
        v_conv = a.v[xCenter, yCenter]
        X, Y, Uw, Vw = tools.window(a,xCenter,yCenter,dist)
        model = fit(model[0], model[1], X, Y, model[2], model[3], Uw, Vw, u_conv, v_conv)
        uMod, vMod = velocity_model(model[0], model[1], model[2], model[3], u_conv, v_conv,X,Y)
        corr = correlation_coef(Uw,Vw,uMod,vMod)
        
        #print('dist:',dist,'Radius',round(model[0],3),'Gamma',
        #      round(model[1],3),'corr',round(corr,3),'x',model[2],
        #      'y',model[3],'u_conv',u_conv,'v_conv',v_conv,
        #      'xC',xCenter,'yC',yCenter) 
    return xCenter, yCenter, model[1], model[0], corr, dist, model[2], model[3], u_conv, v_conv

def fit(coreR, gamma, x, y, fxCenter, fyCenter, Uw, Vw, u_conv, v_conv):
    x = x.ravel()
    y = y.ravel()
    Uw = Uw.ravel()
    Vw = Vw.ravel()
    
    def fun(fitted): #fitted[0]=coreR, fitted[1]=gamma, fitted[2]=fxCenter, fitted[3]=fyCenter, fitted[4]=u_conv, fitted[5]=v_conv
        r = np.hypot(x-fitted[2], y-fitted[3])
        expr2 = np.exp(-r**2/fitted[0]**2)
        z = fitted[1]/(2*np.pi*r) * (1 - expr2)
        z = np.nan_to_num(z)
        zx = (-z + u_conv)*(x-fitted[2]) -Uw
        zy = (z + v_conv)*(y-fitted[3]) -Vw
        zt = np.append(zx,zy)
        return zt
    bnds=([0.001,-100,fxCenter-0.03,fyCenter-0.03],
          [2.0,100,fxCenter+0.03,fyCenter+0.03])
    sol = optimize.least_squares(fun, [coreR,gamma,fxCenter,fyCenter],bounds=bnds,method='dogbox')     
    #Levenberg working!
    #sol = optimize.least_squares(fun, [coreR,gamma,fxCenter,fyCenter],method='lm')
    return sol.x
