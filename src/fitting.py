import numpy as np
from scipy.signal import correlate2d
from scipy import optimize
from scipy.stats import pearsonr

import tools

def correlation_coef(Uw,Vw,u,v):
    #corr_x = (np.mean(np.dot(Uw,u))/(np.sqrt(np.mean(np.dot(Uw,Uw))*np.sqrt(np.mean(np.dot(u,u))))))**0.5
    #corr_y = (np.mean(np.dot(Vw,v))/(np.sqrt(np.mean(np.dot(Vw,Vw))*np.sqrt(np.mean(np.dot(v,v))))))**0.5
    #corr_x = (np.mean(Uw*u)/(np.sqrt(np.mean(Uw**2))*np.sqrt(np.mean(u**2))))**0.5
    #corr_y = (np.mean(Vw*v)/(np.sqrt(np.mean(Vw**2))*np.sqrt(np.mean(v**2))))**0.5
    #R = corr_x*corr_y
    Uw2 = Uw.ravel()
    u2 = u.ravel()
    Vw2 = Vw.ravel()
    v2 = v.ravel()
    R2x = pearsonr(Uw2,u2)
    R2y = pearsonr(Vw2,v2)
    R2 = R2x[0]*R2y[0]
    #print(R2)
    
    return R2

def velocity_model(coreR, gamma, u_conv, v_conv,fxCenter,fyCenter,x,y):
    r = np.hypot(x-fxCenter, y-fyCenter)
    vel = (gamma/(2 * np.pi * r)) * (1 - np.exp(-(r**2)/(coreR)**2))
    vel = np.nan_to_num(vel)
    velx = (vel + u_conv)*(-x+fxCenter)
    vely = (vel + v_conv)*(y-fyCenter)
    return velx, vely

def full_fit(a, xCenter, yCenter, gamma):
    u_conv = a.u[xCenter, yCenter]
    v_conv = a.v[xCenter, yCenter]
    fxCenter = a.dx[xCenter]
    fyCenter = a.dy[yCenter]
    coreR = 0.1
    corrOld = 0.0
    corr = 0.001
    dist = 5

    #print('xC:',xCenter,'yC:',yCenter, 'vort:',gamma)
    while (corr > corrOld):
        corrOld = corr
        coreROld = coreR
        gammaOld = gamma
        dx = a.dx[xCenter+1]-a.dx[xCenter]
        dy = a.dy[yCenter+1]-a.dx[yCenter]
        shiftxCenter = (a.dx[xCenter] -fxCenter)/dx
        shiftyCenter = (a.dy[yCenter] -fyCenter)/dy
        newxCenter = int(shiftxCenter)
        newyCenter = int(shiftyCenter)
        xCenter = xCenter + newxCenter
        yCenter = yCenter + newyCenter
        X, Y, Uw, Vw = tools.window(a,xCenter,yCenter,dist)
        model = [[],[],[],[],[],[]]
        model = fit(X, Y, fxCenter, fyCenter, Uw, Vw, u_conv, v_conv, gamma)
        uMod, vMod = velocity_model(model[0], model[1], model[2], model[3], u_conv, v_conv,X,Y)
        corr = correlation_coef(Uw,Vw,uMod,vMod)
        dist += 1
        #print('dist:',dist-1,'Radius',round(coreR,3),'Gamma',
        #      round(gamma,3),'corr',round(corr,3),'x',fxCenter,
        #      'y',fyCenter,'u_conv',u_conv,'v_conv',v_conv)
           
    return coreROld, gammaOld, corrOld, dist-2, fxCenter, fyCenter, u_conv, v_conv

def fit(x, y, fxCenter, fyCenter, Uw, Vw, u_conv, v_conv, gamma):
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
    bnds=([0.01,-100,fxCenter-0.1,fyCenter-0.1],
          [2.00,+100,fxCenter+0.1,fyCenter+0.1])
    sol = optimize.least_squares(fun, [0.05,gamma,fxCenter,fyCenter],bounds=bnds,method='dogbox')     
    return sol.x
