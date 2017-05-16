import numpy as np

def correlation(a,b):
    print('todo')
    u_i = u_piv - u_c
    u_m = u_model - u_c
    R = np.sqrt( np.mean(((u_i)*(u_m)))/
                ((np.sqrt(np.mean(u_i**2))*np.sqrt(np.mean(u_m**2)))) )
    return R


def model_oseen(a, x, y, coreR):
    r = np.hypot(x, y)
    gamma = 1.0
    u_c = 0.0
    vel = u_c + (gamma/(2 * np.pi * r)) * (1 - np.exp(-(r/coreR)**2))
    return y * vel, -x * vel

