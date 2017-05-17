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
    plt.contourf(field)#, cmap="Greys_r")
    plt.scatter(dirL[1],dirL[0],s=dirL[2],edgecolor='G',facecolor='none')
    plt.scatter(dirR[1],dirR[0],s=dirR[2],edgecolor='Y',facecolor='none')
    #plt.imshow(field, cmap="Greys_r",origin="lower")
    plt.tight_layout()
    
    plt.show()

def plot_quiver(a, xCenter, yCenter, dist, field):
    X, Y = np.meshgrid(a.dx[xCenter-dist:xCenter+dist],
                       a.dy[yCenter-dist:yCenter+dist])
    Uw = a.u[xCenter-dist:xCenter+dist,yCenter-dist:yCenter+dist]
    Vw = a.v[xCenter-dist:xCenter+dist,yCenter-dist:yCenter+dist]

    plt.figure()
    plt.title('Velocity vectors centered at max swirling strength')
    plt.contourf(field[xCenter-dist:xCenter+dist,yCenter-dist:yCenter+dist],
                 extent=[X[0][0], X[0][-1], Y[0][0], Y[-1][0]],origin='lower')
    s = 1
    plt.quiver(X[::s,::s],Y[::s,::s],Vw[::s,::s],Uw[::s,::s])
    
    plt.show()
    
def plot_corr(X, Y, Uw, Vw, uMod, vMod):
    plt.figure()
    plt.title('Correlation')
    s = 1
    plt.quiver(X[::s,::s], Y[::s,::s], Vw[::s,::s],Uw[::s,::s],color='r',scale=15)
    plt.quiver(X[::s,::s], Y[::s,::s], uMod[::s,::s], vMod[::s,::s],color='b',scale=15)
    
    plt.show()
