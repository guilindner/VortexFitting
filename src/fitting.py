import numpy as np

def correlation(Uw,Vw,u,v):
    u = np.nan_to_num(u)
    v = np.nan_to_num(v)
    x_corr = (np.corrcoef(Vw,u))
    y_corr = (np.corrcoef(Uw,v))
    R = (x_corr+y_corr)/2.
    #print('!!!Vw',Uw)
    #print('!!!v',v)
    
    #print('Uw',Vw)
    #print('u',u)
    
    print('x:',x_corr[0])
    print('y:',y_corr[0])
    #sumsquare_z_x = 0.0
    #sumsquare_z_y = 0.0
    #sumsquare_o_x = 0.0
    #sumsquare_o_y = 0.0
    #sum_zo_x = 0.0
    #sum_zo_y = 0.0
    #u = Vw
    #sumsquare_z_x += np.sum(Vw**2)
    #sumsquare_z_y += np.sum(Uw**2)
    #sumsquare_o_x += np.sum(u**2)
    #sumsquare_o_y += np.sum(v**2)
    #sum_zo_x += np.sum(Uw*u)
    #sum_zo_y += np.sum(Vw*v)
    ##sum_zo = sum_zo_x + sum_zo_y
    #corr_x = sum_zo_x**2/(sumsquare_z_x*sumsquare_o_x)
    #corr_y = sum_zo_y**2/(sumsquare_z_y*sumsquare_o_y)  
    #print('corr x',corr_x)
    #print('corr y',corr_y)
    #R = 1.
    
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
    
