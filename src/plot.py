import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def plot_fields(a,vorticity):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)#, sharex='col', sharey='row')
    ax1.imshow(a.u, cmap='seismic')
    ax1.set_title('Velocity u (velocity_s)')
    
    ax2.imshow(a.v, cmap='seismic')
    ax2.set_title('Velocity v (velocity_n)')
    
    
    ax3.imshow(a.w, cmap='seismic')
    ax3.set_title('Velocity w (velocity_z)')
    
    totalvel = np.sqrt(a.u**2 + a.v**2 + a.w**2)
    ax4.set_title('Total velocity (u, v and w)')
    ax4.imshow(vorticity)#, cmap='seismic')
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
    
def plot_corr(a, xCenter, yCenter, dist, gamma):
    X, Y = np.meshgrid(a.dx[xCenter-dist:xCenter+dist],
                       a.dy[yCenter-dist:yCenter+dist])
    Uw = a.u[xCenter-dist:xCenter+dist,yCenter-dist:yCenter+dist]
    Vw = a.v[xCenter-dist:xCenter+dist,yCenter-dist:yCenter+dist]

    plt.figure()
    plt.title('Correlation')
    s = 1
    plt.quiver(X[::s,::s],Y[::s,::s],Vw[::s,::s],Uw[::s,::s],color='r',scale=15)
    #OSEEN
    #gamma = 5.9
    coreR = 0.3
    def velocity(x, y, gamma, coreR):
        r = np.hypot(x-a.dx[xCenter], y-a.dy[yCenter])
        vel = (gamma/(2 * np.pi * r)) * (1 - np.exp(-(r**2)/(coreR)**2))
        return (y-a.dy[yCenter]) * vel, (-x+a.dx[xCenter]) * vel
    
    
    u, v = velocity(X, Y, gamma, coreR)
    plt.quiver(X[::s,::s], Y[::s,::s], u, v,color='b',scale=15)
    u = np.nan_to_num(u)
    v = np.nan_to_num(v)
    #corr = (np.mean(((U)*(u)))/(np.mean(((U**2)**0.5))*np.mean(((u**2)**0.5))))**0.5
    #print(Uw)
    #print(Uw)
    #u = Uw
    #v = Vw
    Uw = Uw - a.u[xCenter,yCenter]
    sumsquare_z_x = 0.0
    sumsquare_z_y = 0.0
    sumsquare_o_x = 0.0
    sumsquare_o_y = 0.0
    sum_zo_x = 0.0
    sum_zo_y = 0.0
    for i in range(len(u[0])):
        for j in range(len(u[1])):
            sumsquare_z_x += (Uw[i,j]-a.u[xCenter,yCenter])**2
            sumsquare_z_y += (Vw[i,j]-a.v[xCenter,yCenter])**2
            sumsquare_o_x += (u[i,j]-a.u[xCenter,yCenter])**2
            sumsquare_o_y += (v[i,j]-a.v[xCenter,yCenter])**2
            sum_zo_x += (Uw[i,j]-a.u[xCenter,yCenter])*(u[i,j]-a.u[xCenter,yCenter])
            sum_zo_y += (Vw[i,j]-a.v[xCenter,yCenter])*(u[i,j]-a.v[xCenter,yCenter])
    sum_zo = sum_zo_x + sum_zo_y
    sumsquare_z = sumsquare_z_x + sumsquare_z_y
    sumsquare_o = sumsquare_o_x + sumsquare_o_y
    correlation = (sum_zo**2)/(sumsquare_z*sumsquare_o)
    print(correlation)
    print('x',np.nanmean(np.corrcoef(Uw.ravel(),u.ravel())))
    print('y',np.nanmean(np.corrcoef(Vw.ravel(),v.ravel())))
    Uw = Uw.ravel()
    #print(Uw)

    
    plt.show()
