import numpy as np
from scipy.signal import correlate2d
from scipy import optimize
from scipy.stats import pearsonr

import tools

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
    
    #print(np.mean(np.dot(Uw,u)))
    #print('correlation x',corr_x)
    #print('correlation y',corr_y)
    #print('corr x',corr_x,'corr y',corr_y,'total correlation',R)
    return R

def velocity_model(a, x, y,xCenter,yCenter, gamma, coreR):
    r = np.hypot(x-a.dx[xCenter], y-a.dy[yCenter])
    vel = (gamma/(2 * np.pi * r)) * (1 - np.exp(-(r**2)/(coreR)**2))
    vel = np.nan_to_num(vel)
    u_conv = a.u[xCenter,yCenter] 
    v_conv = a.v[xCenter,yCenter]
    velx = (vel + u_conv)*(-x+a.dx[xCenter])
    vely = (vel + v_conv)*(y-a.dy[yCenter])
    return velx, vely

def velocity_modelf(a, x, y,xCenter,yCenter,fxCenter,fyCenter, gamma, coreR):
    r = np.hypot(x-fxCenter, y-fyCenter)
    vel = (gamma/(2 * np.pi * r)) * (1 - np.exp(-(r**2)/(coreR)**2))
    vel = np.nan_to_num(vel)
    u_conv = a.u[xCenter,yCenter] 
    v_conv = a.v[xCenter,yCenter]
    velx = (vel + u_conv)*(-x+fxCenter)
    vely = (vel + v_conv)*(y-fyCenter)
    return velx, vely

def full_fit(a, xCenter, yCenter, gamma):
    u_conv = a.u[xCenter, yCenter]
    v_conv = a.v[xCenter, yCenter]
    coreR = 0.1
    corrOld = 0.0
    corr = 0.001
    dist = 3
    #print('xC:',xCenter,'yC:',yCenter, 'vort:',gamma)
    while (corr > corrOld):
        corrOld = corr
        coreROld = coreR
        gammaOld = gamma
        X, Y, Uw, Vw = tools.window(a,xCenter,yCenter,dist)
        fxCenter = a.dx[xCenter]
        fyCenter = a.dy[yCenter]
        #2coreR, gamma = fit(a, X, Y, fxCenter, fyCenter, Uw, Vw, u_conv, v_conv, gamma)
        coreR, gamma, fxCenter, fyCenter = fit4(a, X, Y, fxCenter, fyCenter, Uw, Vw, u_conv, v_conv, gamma)
        uMod, vMod = velocity_model(a, X, Y,xCenter,yCenter, gamma, coreR)
        corr = correlation_coef(Uw,Vw,uMod,vMod)
        dist += 1
        #print('dist:',dist-1,'O Radius',round(coreROld,3),
        #      'N Radius',round(coreR,3),'O Gamma',round(gammaOld,3),
        #      'N Gamma',round(gamma,3),'O corr',round(corrOld,3),
        #      'N corr',round(corr,3))
        
    #2return coreROld, corrOld, dist-2    
    return coreROld, corrOld, dist-2, fxCenter, fyCenter
        
  
def fit(a, x, y, fxCenter, fyCenter, Uw, Vw, u_conv, v_conv, gamma):
    x = x.ravel()
    y = y.ravel()
    Uw = Uw.ravel()
    Vw = Vw.ravel()
    
    def fun(fitted): #fitted[0]=coreR, fitted[1]=gamma
        r = np.hypot(x-fxCenter, y-fyCenter)
        expr2 = np.exp(-r**2/fitted[0]**2)
        z = fitted[1]/(2*np.pi*r) * (1 - expr2)
        z = np.nan_to_num(z)
        zx = (-z + u_conv)*(x-fxCenter) -Uw
        zy = (z + v_conv)*(y-fyCenter) -Vw
        zt = np.append(zx,zy)
        return zt       
    
    sol = optimize.least_squares(fun, [0.5,gamma],method='lm')
    return sol.x
    
def fit4(a, x, y, fxCenter, fyCenter, Uw, Vw, u_conv, v_conv, gamma):
    x = x.ravel()
    y = y.ravel()
    Uw = Uw.ravel()
    Vw = Vw.ravel()
    
    def fun(fitted): #fitted[0]=coreR, fitted[1]=gamma, fitted[2]=xCenter, fitted[3]=yCenter
        r = np.hypot(x-fitted[2], y-fitted[3])
        expr2 = np.exp(-r**2/fitted[0]**2)
        z = fitted[1]/(2*np.pi*r) * (1 - expr2)
        z = np.nan_to_num(z)
        zx = (-z + u_conv)*(x-fitted[2]) -Uw
        zy = (z + v_conv)*(y-fitted[3]) -Vw
        zt = np.append(zx,zy)
        return zt
    bnds=([0.01,-50,fxCenter-0.2,fyCenter-0.2],[2.0,50,
           fxCenter+0.2,fyCenter+0.2])
    sol = optimize.least_squares(fun, [0.5,gamma,fxCenter,fyCenter],bounds=bnds,method='dogbox')
          
    return sol.x
