import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import tools

def plot_fields(a,vorticity):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)#, sharex='col', sharey='row')
    ax1.imshow(a.u, cmap='seismic',origin="lower")
    ax1.set_title('Velocity u (velocity_s)')
    
    ax2.imshow(a.v, cmap='seismic',origin="lower")
    ax2.set_title('Velocity v (velocity_n)')
    
    
    ax3.imshow(a.w, cmap='seismic',origin="lower")
    ax3.set_title('Velocity w (velocity_z)')
    
    totalvel = np.sqrt(a.u**2 + a.v**2 + a.w**2)
    ax4.set_title('Total velocity (u, v and w)')
    ax4.imshow(vorticity,origin="lower")#, cmap='seismic')
    plt.tight_layout()
    
    plt.show()
    
def plot_detection(dirL,dirR,field):
    plt.subplot()
    plt.title('detection')
    plt.contourf(field, cmap="Greys_r")
    plt.scatter(dirL[1],dirL[0],s=dirL[2]*10,edgecolor='G',facecolor='none',label='left')
    plt.scatter(dirR[1],dirR[0],s=dirR[2]*10,edgecolor='Y',facecolor='none',label='right')
    plt.legend()
    #plt.imshow(field, cmap="Greys_r",origin="lower")
    plt.tight_layout()
    
    plt.show()

def plot_quiver(X, Y, Uw, Vw, field):
    plt.figure()
    plt.title('Velocity vectors centered at max swirling strength')
    plt.contourf(field,
                 extent=[X[0][0], X[0][-1], Y[0][0], Y[-1][0]])
    s = 1
    plt.quiver(X[::s,::s],Y[::s,::s],Vw[::s,::s],Uw[::s,::s])
    
    plt.show()
    
def plot_corr(X, Y, Uw, Vw, uMod, vMod):
    plt.figure()
    plt.title('Correlation')
    s = 1
    if (X.size > 400):
        s = 2
    plt.quiver(X[::s,::s], Y[::s,::s], Vw[::s,::s],Uw[::s,::s],
               color='r',label='data')
    plt.quiver(X[::s,::s], Y[::s,::s], vMod[::s,::s], uMod[::s,::s],
               color='b',label='model')
    plt.legend()
    plt.show()
