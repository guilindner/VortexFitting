import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import norm
import re

import tools
import fitting

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

def plot_fit(X, Y, Uw, Vw, uMod, vMod, xc, yc, coreR, gamma, u_conv, v_conv, corr,i,j):
    plt.figure()
    s = 1
    if (X.size > 400):
        s = 1
    plt.quiver(X[::s,::s], Y[::s,::s], Uw[::s,::s],Vw[::s,::s],
               color='r',label='data')
    plt.quiver(X[::s,::s], Y[::s,::s], uMod[::s,::s], vMod[::s,::s],
               color='b',label='model', alpha=0.5)
    circle1=plt.Circle((xc,yc),coreR,color='k',alpha=0.05)
    plt.gca().add_artist(circle1)
    plt.gca().scatter([xc], [yc], marker='+', color='k', s=100)
    plt.legend()
    plt.grid()
    plt.axes().set_aspect('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(r'r=%s $\Gamma$=%s u=%s v=%s C=%s' %(round(coreR,2),round(gamma,2),round(u_conv,2),round(v_conv,2),round(corr,2)))
    plt.savefig('../results/vortex%i_%i.png' %(i,j),format='png')
    plt.close('all')

def plot_fit_test(X, Y, Uw, Vw, uMod, vMod, xc, yc, coreR, gamma, u_conv, v_conv, corr):
    plt.figure()
    s = 1
    if (X.size > 400):
        s = 1
    plt.quiver(X[::s,::s], Y[::s,::s], Uw[::s,::s],Vw[::s,::s],
               color='r',label='data')
    plt.quiver(X[::s,::s], Y[::s,::s], uMod[::s,::s], vMod[::s,::s],
               color='b',label='model', alpha=0.5)
    circle1=plt.Circle((xc,yc),coreR,color='k',alpha=0.05)
    plt.gca().add_artist(circle1)
    plt.gca().scatter([xc], [yc], marker='+', color='k', s=100)
    plt.legend()
    plt.grid()
    plt.axes().set_aspect('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(r'r=%s $\Gamma$=%s u=%s v=%s C=%s' %(round(coreR,2),round(gamma,2),round(u_conv,2),round(v_conv,2),round(corr,2)))
    plt.show()

def plot_accepted(a,vortices,field):
    plt.subplot()
    plt.imshow(field, origin='lower', cmap="Greys_r")
    plt.xlabel('x')
    plt.ylabel('y')
    dx = a.dx[5]-a.dx[4]
    dy = a.dy[5]-a.dy[4]
    for i,line in enumerate(vortices):
        if vortices[i][1] > 0:
            orient = 'Y'
        else:
            orient = 'Y'
        circle1=plt.Circle((line[2]/dx,line[3]/dy),line[0]/dx,
                            edgecolor=orient,facecolor='none',gid='vortex%i' % i)
        plt.gca().add_artist(circle1)

    ##Comparing data
    #fileIn = open('../data/dazin.dat', 'r')
    #for line in fileIn:
        #xComp = int(float(line.split()[1]))
        #yComp = int(float(line.split()[2]))
        #gammaComp = float(line.split()[3])
        #rComp = float(line.split()[4])
        #if gammaComp > 0:
            #orient = 'R'
        #else:
            #orient = 'R'
        #circle2=plt.Circle((xComp,yComp),rComp,edgecolor=orient,facecolor='none')
        #plt.gca().add_artist(circle2)

    #plt.legend()
    plt.tight_layout()
    plt.savefig('../results/accepted.svg', format='svg')
    plt.savefig('../results/tk.png', format='png', transparent=True)
    create_links('../results/accepted.svg',vortices)
    #plt.show()

def plot_vortex(a,vortices):
    outfile = open('../results/vortices.dat','w')
    outfile.write('radius gamma X Y u_c v_c dist corr\n')
    for i,line in enumerate(vortices):
        #print(line)
        outfile.write("{0} {1} {2} {3} {4} {5} {6} {7} \n".format(line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7]))
        print('r:',line[0],'gamma:',line[1], 'x:',line[2],
         'y',line[3],'dist',line[6],'corr',line[7])
        dx = a.dx[5]-a.dx[4]
        dy = a.dy[5]-a.dy[4]
        X, Y, Uw, Vw = tools.window(a,round(line[2]/dx,0),round(line[3]/dy,0),line[6])
        uMod, vMod = fitting.velocity_model(line[0], line[1],
         line[2], line[3], line[4], line[5], X, Y)
        corr = fitting.correlation_coef(Uw,Vw,uMod,vMod)
        plot_fit(X, Y, Uw, Vw, uMod, vMod, line[2],line[3], line[0], line[1], line[4], line[5], corr,i,1)
        corr = fitting.correlation_coef(Uw-line[4],Vw-line[5],uMod-line[4],vMod-line[5])
        plot_fit(X, Y, Uw-line[4], Vw-line[5], uMod-line[4], vMod-line[5], line[2],
                   line[3], line[0], line[1], line[4], line[5], corr,i,2)

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
            fileOut.write('   <a href="vortex%i_1.png">\n' % i)
            fileOut.write(line)
            fileOut.write('   <title>Vortex %i: r = %s gamma = %s</title>\n' % (i,round(vortices[i][3],1),round(vortices[i][2],1)) )
            i = i + 1
            vortex_found = True
        else:
            fileOut.write(line)


