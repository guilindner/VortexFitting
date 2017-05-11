import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def plot_fields(a):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)#, sharex='col', sharey='row')
    ax1.imshow(a.u, cmap='seismic')
    ax1.set_title('Velocity U (velocity_s)')
    
    ax2.imshow(a.v, cmap='seismic')
    ax2.set_title('Velocity V (velocity_n)')
    
    
    ax3.imshow(a.w, cmap='seismic')
    ax3.set_title('Velocity W (velocity_z)')
    
    totalvel = np.sqrt(a.v**2 + a.w**2)
    ax4.set_title('Total velocity (v and w)')
    ax4.imshow(totalvel, cmap='seismic')
    plt.tight_layout()
    
    plt.show()
    
def plot_detection(dirL,dirR,field):
    plt.subplot()
    plt.title('detection')
    plt.contourf(field, cmap="Greys_r")
    plt.scatter(dirL[1],dirL[0],s=dirL[2]*10,edgecolor='G',facecolor='none')
    plt.scatter(dirR[1],dirR[0],s=dirR[2]*10,edgecolor='Y',facecolor='none')
    
    #plt.imshow(field, cmap="Greys_r",origin="lower")
    plt.tight_layout()
    
    plt.show()

def plot_quiver(a, xCenter, yCenter, dist, field):
    X, Y = np.meshgrid(a.dx[xCenter-dist:xCenter+dist],
                       a.dy[yCenter-dist:yCenter+dist])
    U = a.u[xCenter-dist:xCenter+dist,yCenter-dist:yCenter+dist]
    V = a.v[xCenter-dist:xCenter+dist,yCenter-dist:yCenter+dist]
    plt.figure()
    plt.title('Velocity vectors centered at max swirling strength')
    plt.contourf(field[xCenter-dist:xCenter+dist,yCenter-dist:yCenter+dist],
                 extent=[X[0][0], X[0][-1], Y[0][0], Y[-1][0]])
    s = 2
    plt.quiver(X[::s,::s], Y[::s,::s], U[::s,::s], V[::s,::s])
    
    plt.show()
