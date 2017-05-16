import numpy as np
from scipy import optimize
np.seterr(divide='ignore', invalid='ignore') #when r = to 0, dont give error

#todo oseen_vorticity?

def velocity(x, y, coreR):
    r = np.hypot(x, y)
    gamma = 1.0
    vel = (gamma/(2 * np.pi * r)) * (1 - np.exp(-(r**2)/(coreR)**2))
    return y * vel, -x * vel


def main():
    import matplotlib.pyplot as plt

    dom = np.linspace(-0.02, 0.02, 16)
    x, y = np.meshgrid(dom, dom)
    x_vel = x[::1, ::1] 
    y_vel = y[::1, ::1]
    
    coreR = 1.0
    u, v = velocity(x_vel, y_vel, coreR)
    plt.quiver(x_vel, y_vel, u, v)


    local_vorticity = 1.0
    
    def jac_oseen(radius):
        return 1.0 
    def oseenx(radius):
        x = local_vorticity/(2*np.pi*(radius-velx))*(1-np.exp(-((radius-velx)/coreR)**2))
        #y = local_vorticity/(2*np.pi*(radius-vely))*(1-np.exp(-((radius-vely)/coreR)**2))
        return x
    def oseeny(radius):
        y = local_vorticity/(2*np.pi*(radius-vely))*(1-np.exp(-((radius-vely)/coreR)**2))
        return y
    for i in range(2):
        for j in range(2):
            velx = u[i][j]
            vely = v[i][j]
            print(u[i][j],v[i][j])
            root_x = optimize.root(oseenx,0.0,method='lm')
            root_y = optimize.root(oseeny,0.0,method='lm')
            print(root_x.x, root_y.x)

    
    plt.show()


if __name__ == '__main__':
    main()
