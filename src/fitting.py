import numpy as np
from scipy.signal import correlate2d
from scipy import optimize

import tools

def correlation_coef(Uw,Vw,u,v):
    corr_x = np.mean(Uw*u)/(np.sqrt(np.mean(Uw**2))*np.sqrt(np.mean(u**2)))
    corr_y = np.mean(Vw*v)/(np.sqrt(np.mean(Vw**2))*np.sqrt(np.mean(v**2)))  
    #print('corr x',corr_x)
    #print('corr y',corr_y)
    R = corr_x*corr_y
    return R

def model_oseen_x(coreR):
    #xCenter = 37
    #yCenter = 47
    #gamma = 22.4031
    #x, y, Uw, Vw = tools.window(a,xCenter,yCenter,dist)
    #r = np.hypot(x-a.dx[xCenter], y-a.dy[yCenter])
    r = np.hypot(2-a.dx[xCenter], 2-a.dy[yCenter])
    velx = a.u[x,y] - (gamma/(2 * np.pi * r)) * (1 - np.exp(-(r**2)/(coreR)**2))
    vel = np.nan_to_num(vel)
    u_conv = a.u[xCenter,yCenter] 
    velx2 = (velx + u_conv)*(-x+a.dx[xCenter])
    return velx2

def velocity_model(a, x, y,xCenter,yCenter, gamma, coreR):
    r = np.hypot(x-a.dx[xCenter], y-a.dy[yCenter])
    vel = (gamma/(2 * np.pi * r)) * (1 - np.exp(-(r**2)/(coreR)**2))
    vel = np.nan_to_num(vel)
    u_conv = a.u[xCenter,yCenter] 
    v_conv = a.v[xCenter,yCenter]
    velx = (vel + u_conv)*(-x+a.dx[xCenter])
    vely = (vel + v_conv)*(y-a.dy[yCenter])
    return velx, vely
    
    
 
def super_fitx(a, x, y, xCenter, yCenter, Uw, velx, gamma):
    x = x.ravel()
    y = y.ravel()
    Uw = Uw.ravel()
    def funx(coreR):
        r = np.hypot(x-a.dx[xCenter], y-a.dy[yCenter])
        expr2 = np.exp(-r**2/coreR**2)
        z = -gamma/(2*np.pi*r) * (1 - expr2)
        z = np.nan_to_num(z)
        z = (z + velx)*(-x+a.dx[xCenter]) +Uw
        #print(Uw.reshape(10,10))
        return z        
  
    sol = optimize.root(funx, 0.1, jac=False, method='lm')
    #print(sol.fun.reshape(10,10))
    return sol.x

def super_fity(a, x, y, xCenter, yCenter, Vw, vely, gamma):
    x = x.ravel()
    y = y.ravel()
    Vw = Vw.ravel()
    def funy(coreR):
        r = np.hypot(x-a.dx[xCenter], y-a.dy[yCenter])
        expr2 = np.exp(-r**2/coreR**2)
        z = gamma/(2*np.pi*r) * (1 - expr2)
        z = np.nan_to_num(z)
        z = (z + vely)*(y-a.dy[yCenter]) -Vw
        return z         

    sol = optimize.root(funy, 0.1, jac=False, method='lm')
    return sol.x
