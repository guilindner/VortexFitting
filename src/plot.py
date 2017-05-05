import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def plot_fields(a):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)#, sharex='col', sharey='row')
    ax1.imshow(a.u, cmap='seismic')
    ax1.set_title('Velocity U (velocity_s)')
    
    ax2.imshow(a.v, cmap='seismic')
    ax2.set_title('Velocity V (velocity_n)')
    
    totalvel = np.sqrt(a.u**2 + a.v**2)
    ax3.set_title('Total velocity')
    ax3.imshow(totalvel, cmap='seismic')
    
    ax4.set_title('Total velocity')
    ax4.imshow(totalvel, cmap='seismic')
    plt.tight_layout()
    
    plt.show()
    
def plot_detection(dirL,dirR,field):
    plt.subplot()
    plt.title('Vorticity, rotation')
    plt.scatter(dirL[0],dirL[1],s=dirL[2],edgecolor='G',facecolor='none')
    plt.scatter(dirR[0],dirR[1],s=dirR[2],edgecolor='Y',facecolor='none')
    plt.imshow(field, vmax=10, cmap="Greys_r")
    plt.tight_layout()
    
    plt.show()

def plot_quiver(a, xCenter, yCenter, dist=15):
    X, Y = np.meshgrid(a.dx[xCenter-dist:xCenter+dist],
                       a.dy[yCenter-dist:yCenter+dist])
    U = a.u[yCenter-dist:yCenter+dist,xCenter-dist:xCenter+dist]
    V = a.v[yCenter-dist:yCenter+dist,xCenter-dist:xCenter+dist]
    
    plt.figure()
    plt.title('Arrows scale with plot width, not view')
    Q = plt.quiver(X, -Y, U, V)
 
    plt.show()
