import numpy as np
from scipy.signal import correlate2d

def correlation(Uw,Vw,u,v):
    #Uw[0] = (Uw[0]-np.mean(Uw[0]))/np.std(Uw[0])
    #u[0] = (u[0]-np.mean(u[0]))/np.std(u[0])
    #x_corr = (correlate2d(Uw,u))
    #y_corr = (correlate2d(Vw,v))
    #R = (x_corr+y_corr)/2.
    #print('!!!Vw',Uw)
    #print('!!!v',v)
    
    #print('Uw',Vw)
    #print('u',u)
    
    
    corr_x = np.mean(Uw*u)/(np.sqrt(np.mean(Uw**2))*np.sqrt(np.mean(u**2)))
    corr_y = np.mean(Vw*v)/(np.sqrt(np.mean(Vw**2))*np.sqrt(np.mean(v**2)))  
    print('corr x',corr_x)
    print('corr y',corr_y)
    R = (corr_x+corr_y)/2
    print(R)
    return R

def Ruim_model_oseen(a, x, y, coreR):
    r = np.hypot(x, y)
    gamma = 1.0
    u_c = 0.0
    vel = u_c + (gamma/(2 * np.pi * r)) * (1 - np.exp(-(r/coreR)**2))
    return y * vel, -x * vel

def velocity_model(a, x, y,xCenter,yCenter, gamma, coreR):
    r = np.hypot(x-a.dx[xCenter], y-a.dy[yCenter])
    vel = (gamma/(2 * np.pi * r)) * (1 - np.exp(-(r**2)/(coreR)**2))
    vel = np.nan_to_num(vel)
    u_conv = a.u[xCenter,yCenter] 
    v_conv = a.v[xCenter,yCenter]
    velx = (vel + u_conv)*(-x+a.dx[xCenter])
    vely = (vel + v_conv)*(y-a.dy[yCenter])
    return velx, vely
    
