import numpy as np
from scipy.signal import correlate2d
from scipy import optimize
from scipy.stats import pearsonr

import tools
import plot

def correlation_coef(Uw,Vw,u,v):
    #corr_x = (np.mean(np.dot(Uw,u))/(np.sqrt(np.mean(np.dot(Uw,Uw))*np.sqrt(np.mean(np.dot(u,u))))))**0.5
    #corr_y = (np.mean(np.dot(Vw,v))/(np.sqrt(np.mean(np.dot(Vw,Vw))*np.sqrt(np.mean(np.dot(v,v))))))**0.5
    corr_x = (np.mean(Uw*u)/(np.sqrt(np.mean(Uw**2))*np.sqrt(np.mean(u**2))))**0.5
    corr_y = (np.mean(Vw*v)/(np.sqrt(np.mean(Vw**2))*np.sqrt(np.mean(v**2))))**0.5
    R = corr_x*corr_y
    Uw2 = Uw.ravel()
    u2 = u.ravel()
    Vw2 = Vw.ravel()
    v2 = v.ravel()
    R2x = pearsonr(Uw2,u2)
    R2y = pearsonr(Vw2,v2)
    R2 = R2x[0]*R2y[0]
    #print(R2)
    
    return R

def velocity_model(coreR, gamma, fxCenter,fyCenter, u_conv, v_conv,x,y):
    r = np.hypot(x-fxCenter, y-fyCenter)
    vel = (gamma/(2 * np.pi * r)) * (1 - np.exp(-(r**2)/(coreR)**2))
    vel = np.nan_to_num(vel)
    velx = (vel + u_conv)*(-x+fxCenter)
    vely = (vel + v_conv)*(y-fyCenter)
    return velx, vely

def full_fit(coreR, gamma, a, xCenter, yCenter):
    model = [[],[],[],[],[],[]]
    model[1] = gamma
    fxCenter = a.dx[xCenter]
    fyCenter = a.dy[yCenter]
    model[2] = fxCenter
    model[3] = fyCenter
    dx = a.dx[xCenter+1]-a.dx[xCenter]
    dy = a.dy[yCenter+1]-a.dx[yCenter]
    model[0] = 0.05
    corrOld = 0.0
    corr = 0.001
    dist = 3
    model[2] = fxCenter
    model[3] = fyCenter
    for i in range(10):
        #print('iter',i)
    #while (corr > corrOld):
        #yCenter = yCenter +1
        #print('x diff', model[2]- fxCenter)
        if (model[2]-fxCenter > dx/2):
            #print('reduce x!')
            xCenter = xCenter -1
        elif (model[2]-fxCenter < -dx/2):
            #print('increase x')
            xCenter = xCenter +1
        fxCenter = model[2]
        #print('y diff',model[3]- fyCenter)
        if (model[3]-fyCenter > dy/2):
            #print('reduce y!')
            yCenter = yCenter -1
        elif (model[3]-fyCenter < -dy/2):
            #print('increase y!')
            yCenter = yCenter +1   
        fxCenterOld = model[2]
        fyCenterOld = model[3]         
        dist = i + 2
        distOld = dist
        corrOld = corr
        coreROld = model[0]
        gammaOld = model[1]
        u_conv = a.u[xCenter, yCenter]
        v_conv = a.v[xCenter, yCenter]
        u_convOld = u_conv
        v_convOld = v_conv
        #print(xCenter, yCenter)
        X, Y, Uw, Vw = tools.window(a,xCenter,yCenter,dist)
        model = fit(model[0], model[1], X, Y, model[2], model[3], Uw, Vw, u_conv, v_conv)
        uMod, vMod = velocity_model(model[0], model[1], model[2], model[3], u_conv, v_conv,X,Y)
        corr = correlation_coef(Uw,Vw,uMod,vMod)
        #print('dist:',dist,'Radius',round(model[0],3),'Gamma',
        #      round(model[1],3),'corr',round(corr,3),'x',model[2],
        #      'y',model[3],'u_conv',u_conv,'v_conv',v_conv)
        print(corr)
        #plot.plot_corr(X, Y, Uw, Vw, uMod, vMod, model[0], corr)
           
    return coreROld, gammaOld, corrOld, distOld, fxCenterOld, fyCenterOld, u_convOld, v_convOld, distOld-1

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
    bnds=([0.01,-100,fxCenter-0.05,fyCenter-0.05],
          [2.00,+100,fxCenter+0.05,fyCenter+0.05])
    sol = optimize.least_squares(fun, [coreR,gamma,fxCenter,fyCenter],bounds=bnds,method='dogbox')     
    return sol.x
