import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import norm

import tools

def plot_fields(a,field):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)#, sharex='col', sharey='row')
    ax1.imshow(a.u, cmap='seismic',origin="lower")
    ax1.set_title('Velocity u (velocity_s)')
    
    ax2.imshow(a.v, cmap='seismic',origin="lower")
    ax2.set_title('Velocity v (velocity_n)')
    
    
    ax3.imshow(a.w, cmap='seismic',origin="lower")
    ax3.set_title('Velocity w (velocity_z)')
    
    totalvel = np.sqrt(a.u**2 + a.v**2 + a.w**2)
    ax4.set_title('Total velocity (u, v and w)')
    ax4.imshow(field,origin="lower", cmap='seismic')
    plt.tight_layout()
    
    plt.show()
    
def plot_fields2(a,vorticity, *args):
    plt.subplot()
    if (args[0] == True):
        a.u = a.u.T
        a.v = a.v.T
        a.w = a.w.T
        vorticity = vorticity.T
    plt.imshow(a.u, cmap='seismic',origin="lower")
    #plt.title('Velocity u')
    plt.show()
    
    plt.imshow(a.v, cmap='seismic',origin="lower")
    #plt.title('Velocity v')
    plt.show()
    
    plt.imshow(a.w, cmap='seismic',origin="lower")
    #plt.title('Velocity w')
    plt.show()
    
    totalvel = np.sqrt(a.u**2 + a.v**2 + a.w**2)
    #plt.title('Total velocity (u, v and w)')
    plt.imshow(vorticity,origin="lower", cmap='seismic')
    plt.show()
    
def plot_detection(dirL,dirR,field, *args):
    plt.subplot()
    if (args[0] == True):
        field = field.T
        plt.scatter(dirL[0],dirL[1],s=dirL[2]*10,edgecolor='G',facecolor='none',label='left')
        plt.scatter(dirR[0],dirR[1],s=dirR[2]*10,edgecolor='Y',facecolor='none',label='right')
    else:
        plt.scatter(dirL[1],dirL[0],s=dirL[2]*10,edgecolor='G',facecolor='none',label='left')
        plt.scatter(dirR[1],dirR[0],s=dirR[2]*10,edgecolor='Y',facecolor='none',label='right')
    
    #plt.title('detection')
    #plt.contourf(field, cmap="Greys_r")

    plt.imshow(field, origin='lower', cmap="Greys_r")
    plt.xlabel('x')
    plt.ylabel('y')
    #plt.legend()
    #plt.imshow(field, cmap="Greys_r",origin="lower")
    plt.tight_layout()
    
    plt.show()

def plot_quiver(X, Y, Uw, Vw, field):
    plt.figure()
    #plt.title('Velocity vectors centered at max swirling strength')
    plt.contourf(field,
                 extent=[X[0][0], X[0][-1], Y[0][0], Y[-1][0]])
    s = 1
    plt.quiver(X[::s,::s],Y[::s,::s],Vw[::s,::s],Uw[::s,::s])
    
    plt.show()
    
def plot_corr(X, Y, Uw, Vw, uMod, vMod, coreR, corr):
    plt.figure()
    plt.title('Correlation')
    s = 1
    if (X.size > 400):
        s = 2
    plt.quiver(X[::s,::s], Y[::s,::s], Uw[::s,::s],Vw[::s,::s],
               color='r',label='data')
    plt.quiver(X[::s,::s], Y[::s,::s], uMod[::s,::s], vMod[::s,::s],
               color='b',label='model')
    plt.legend()
    plt.title('core Radius = %s Correlation = %s' %(round(coreR,3),round(corr,3)))
    plt.show()

def plot_debug(X, Y, Uw, Vw, uMod, vMod, coreR, corr):
    plt.figure()
    plt.title('Correlation')
    s = 1
    if (X.size > 400):
        s = 2
    plt.quiver(X[::s,::s], Y[::s,::s], Uw[::s,::s],Vw[::s,::s],
               color='r',label='data',scale=50)
    plt.quiver(X[::s,::s], Y[::s,::s], uMod[::s,::s], vMod[::s,::s],
               color='b',label='model',scale=50)
    plt.legend()
    plt.show()

def plot_radius(vortices):
    fig, ax = plt.subplots()
    data = []
    for i in range(len(vortices)):
        print(vortices[i][3])
        data.append(vortices[i][3])
    data = np.array(data)
    (mu, sigma) = norm.fit(data)
    print(mu, sigma)
    x = np.linspace(-sigma*3,sigma*3,data.size)
    #x = np.linspace(min(data), max(data), 50)
    y = norm.pdf(x,0.0,sigma)
    ax.plot(x, y, 'r--', linewidth=2)
    ax.set_xlim(-0.3,0.3)
    
    
    #plt.figure()
    #data = []
    #for i in range(len(vortices)):
    #    print(vortices[i][3])
    #    data.append(vortices[i][3])
    #data = np.array(data)
    #print(data)
    #bins = np.arange(min(data), max(data), 0.01)
    #print(bins)
    #plt.hist(data, bins=bins, alpha=0.5)
    
    plt.show()

