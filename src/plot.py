import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import norm
import re

import tools

def plot_fields(a,field):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)#, sharex='col', sharey='row')
    ax1.imshow(a.u, cmap='seismic',origin="lower")
    ax1.set_title('Velocity u (velocity_s)')
    
    ax2.imshow(a.v, cmap='seismic',origin="lower")
    ax2.set_title('Velocity v (velocity_n)')
    
    
    ax3.imshow(a.w, cmap='seismic',origin="lower")
    ax3.set_title('Velocity w (velocity_z)')
    
    ax4.set_title('Vorticity')
    ax4.imshow(field,origin="lower", cmap='seismic')
    plt.tight_layout()
    
    plt.show()
    
def plot_detect(dirL,dirR,field, *args):
    plt.subplot()
    if (args[0] == True):
        field = field.T
        plt.scatter(dirL[0],dirL[1],edgecolor='G',facecolor='G',label='left')
        plt.scatter(dirR[0],dirR[1],edgecolor='Y',facecolor='Y',label='right')
    else:
        plt.scatter(dirL[1],dirL[0],edgecolor='G',facecolor='G',label='left')
        plt.scatter(dirR[1],dirR[0],edgecolor='Y',facecolor='Y',label='right')
    
    plt.title('Detected possible vortices')
    #plt.contourf(field, cmap="Greys_r")

    plt.imshow(field, origin='lower', cmap="Greys_r")
    plt.xlabel('x')
    plt.ylabel('y')
    #plt.legend()
    #plt.imshow(field, cmap="Greys_r",origin="lower")
    plt.tight_layout()
    
    plt.show()

def plot_quiver(X,Y,Uw,Vw,field):
    plt.figure()
    #plt.title('Velocity vectors centered at max swirling strength')
    plt.contourf(field,
                 extent=[X[0][0], X[0][-1], Y[0][0], Y[-1][0]])
    s = 1
    plt.quiver(X[::s,::s],Y[::s,::s],Uw[::s,::s],Vw[::s,::s])
    plt.show()
    
def plot_fit(X, Y, Uw, Vw, uMod, vMod, xc, yc, coreR, gamma, corr,i):
    plt.figure()
    s = 1
    if (X.size > 400):
        s = 1
    plt.quiver(X[::s,::s], Y[::s,::s], Uw[::s,::s],Vw[::s,::s],
               color='r',label='data')
    plt.quiver(X[::s,::s], Y[::s,::s], uMod[::s,::s], vMod[::s,::s],
               color='b',label='model', alpha=0.5)
    circle1=plt.Circle((xc,yc),coreR,color='r',alpha=0.05)
    plt.gca().add_artist(circle1)
    plt.legend()
    plt.grid()
    plt.axes().set_aspect('equal')
    plt.title('Radius = %s Gamma = %s Corr = %s' %(round(coreR,3),round(gamma,3),round(corr,3)))
    plt.savefig('../results/vortex%i.png' % i,format='png')
    plt.close('all')
    
def plot_accepted(vortices,field):
    plt.subplot()
    plt.imshow(field, origin='lower', cmap="Greys_r")
    plt.xlabel('x')
    plt.ylabel('y')
    for i in range(len(vortices)):
        if vortices[i][2] > 0:
            orient = 'Y'
        else:
            orient = 'Y'
        circle1=plt.Circle((vortices[i][0],vortices[i][1]),
                            vortices[i][3],edgecolor=orient,facecolor='none',gid='vortex%i' % i)
        plt.gca().add_artist(circle1)
        
    fileIn = open('../data/dazin.dat', 'r')
    
    for line in fileIn:
        xComp = int(float(line.split()[1]))
        yComp = int(float(line.split()[2]))
        gammaComp = float(line.split()[3])
        rComp = float(line.split()[4])
        if gammaComp > 0:
            orient = 'R'
        else:
            orient = 'R'
        circle2=plt.Circle((xComp,yComp),rComp,edgecolor=orient,facecolor='none')
        plt.gca().add_artist(circle2)       
                
    #plt.legend()
    plt.tight_layout()
    plt.savefig('../results/accepted.svg', format='svg')
    create_links('../results/accepted.svg',vortices)
    #plt.show()
    
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
    
def create_links(path,vortices):
    fileIn = open("../results/accepted.svg","r")
    fileOut = open("../results/linked.svg","w")
    i = 0
    vortex_found = False
    for line in fileIn:
        if "</g>" in line:
            if vortex_found == True:
                fileOut.write(line)
                fileOut.write('   </a>\n')
                vortex_found = False
            else:
                fileOut.write(line)
        elif "vortex" in line:
            fileOut.write('   <a href="vortex%i.png">\n' % i)
            fileOut.write(line)
            fileOut.write('   <title>Vortex %i: r = %s gamma = %s</title>\n' % (i,round(vortices[i][3],3),round(vortices[i][2],3)) )
            i = i + 1
            vortex_found = True
        else:
            fileOut.write(line)
         
    
