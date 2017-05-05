import numpy as np
np.seterr(divide='ignore', invalid='ignore') #when r = to 0, dont give error

#todo oseen_vorticity?

def velocity(x, y, coreR):
    r = np.hypot(x, y)
    gamma = 1.0
    vel = (gamma/(2 * np.pi * r)) * (1 - np.exp(-(r**2)/(coreR)**2))
    return y * vel, -x * vel


def main():
    import matplotlib.pyplot as plt

    dom = np.linspace(-1.0, 1.0, 145)
    x, y = np.meshgrid(dom, dom)
    x_vel = x[::12, ::12] 
    y_vel = y[::12, ::12]
    t = 1.0
    
    coreR = 0.1
    u, v = velocity(x_vel, y_vel, coreR)
    plt.quiver(x_vel, y_vel, u, v)

    plt.show()


if __name__ == '__main__':
    main()
