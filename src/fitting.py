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
    
    
def Bfunx(coreR):
    r2 = (x-xCenter)**2 + (y-yCenter)**2
    expr2 = np.exp(-r2/coreR[0])
    return velx - gamma*(y - yCenter)/r2 * (1 - expr2)

def Bfuny(coreR):
    r2 = (x-xCenter)**2 + (y-yCenter)**2
    expr2 = np.exp(-r2/coreR[0])
    return vely + gamma*(x - xCenter)/r2 * (1 - expr2)

def Bjacx(coreR):
    r2 = (x-xCenter)**2 + (y-yCenter)**2
    expr2 = np.exp(-r2/coreR[0])
    velx_x = -gamma*(y-yCenter)*2.*(x-xCenter)/r2 * (1./r2-(1./r2+1./coreR[0])*expr2)
    velx_y = -gamma*(((y-yCenter)**2-(x-xCenter)**2)/r2**2 * (1.-expr2)-2.*(y-yCenter)**2/(coreR[0]*r2)*expr2)
    return velx_y, velx_x
    
def Bjacy(coreR):
    r2 = (x-xCenter)**2 + (y-yCenter)**2
    expr2 = np.exp(-r2/coreR)
    vely_x = gamma*(((y-yCenter)**2-(x-xCenter)^2)/r2**2 * (1.-expr2)-2.*(y-yCenter)**2/(coreR*r2)*expr2)
    vely_y = gamma*(y-yCenter)*2.*(x-xCenter)/r2 * (1./r2-(1./r2+1./coreR)*expr2)
    return vely_x, vely_y
    
def super_fitx(x,y,velx,vely,gamma):
    x = x.ravel()
    y = y.ravel()
    velx = velx.ravel()
    vely = vely.ravel()
    def funx(coreR):
        r2 = x**2 + y**2
        expr2 = np.exp(-r2/coreR[0])
        return velx - gamma*(y)/r2 * (1 - expr2)
    def jacx(coreR):
        r2 = (x)**2 + (y)**2
        expr2 = np.exp(-r2/coreR[0])
        velx_x = -gamma*(y)*2.*(x)/r2 * (1./r2-(1./r2+1./coreR[0])*expr2)
        velx_y = -gamma*(((y)**2-(x)**2)/r2**2 * (1.-expr2)-2.*(y)**2/(coreR[0]*r2)*expr2)
        return velx_y, velx_x    
    sol = optimize.root(funx, 0.3, jac=jacx, method='lm')
    #print('x,y,vx,vy,gamma',x,y,velx,vely,gamma)
    #print(sol)
    return sol.x

def super_fity(x,y,velx,vely,gamma):
    x = x.ravel()
    y = y.ravel()
    velx = velx.ravel()
    vely = vely.ravel()
    def funy(coreR):
        r2 = x**2 + y**2
        expr2 = np.exp(-r2/coreR[0])
        return vely + gamma*(x)/r2 * (1 - expr2)
    def jacy(coreR):
        r2 = (x)**2 + (y)**2
        expr2 = np.exp(-r2/coreR[0])
        vely_x = gamma*(((y)**2-(x)**2)/r2**2 * (1.-expr2)-2.*(y)**2/(coreR*r2)*expr2)
        vely_y = gamma*(y)*2.*(x)/r2 * (1./r2-(1./r2+1./coreR)*expr2)
        return vely_x, vely_y 
    sol = optimize.root(funy, 0.3, jac=jacy, method='lm')
    #print('x,y,vx,vy,gamma',x,y,velx,vely,gamma)
    return sol.x    
    
