import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy import ndimage
from photutils import find_peaks

#---- LOAD DATA ----#
file1 = 'Data/test_data.nc'
#file2 = 'Data/plan_normal_+1.203_xyz_00.gradient.nc'
#file3 = 'Data/plan_normal_+1.203_xyz.norm.nc'
#file4 = 'Data/plan_normal_+1.203_xyz_00.detection_tourb_adapt.nc'
g1=Dataset(file1,'r')
gv=g1.variables
gd=g1.dimensions
#U = gv['U']
#timestep = gv['time'].size -1
#z,y,x
#print(gv['U'][1,1,100])
#vel = np.array(gv['U'][21,:,:])
velU = np.array(gv['velocity_s'][0])
velV = np.array(gv['velocity_n'][0])
velV = np.einsum('ij->ji',velV)
velU = np.einsum('ij->ji',velU)
print('Shape of vel:',velV.shape)
#x = np.array(vel[])
#sizex = len(gd['resolution_x'])
#sizey = len(gd['resolution_y'])
#sizez = len(gd['resolution_z'])
sizex = len(gd['resolution_n'])
sizey = len(gd['resolution_z'])

#print('Resolution x,y,z:',sizex,sizey,sizez)
print('Resolution x,y',sizex,sizey)

#---- Second-order difference approximation ----#
dx = np.array(gv['grid_n'])
dy = np.array(gv['grid_z'])

dudx = np.zeros((sizex,sizey))
dudy = np.zeros((sizex,sizey))
dudz = np.zeros((sizex,sizey))
dvdx = np.zeros((sizex,sizey))
dvdy = np.zeros((sizex,sizey))
dvdz = np.zeros((sizex,sizey))
dwdx = np.zeros((sizex,sizey))
dwdy = np.zeros((sizex,sizey))
L = np.zeros((sizex,sizey))
ls2 = np.zeros((sizex,sizey))
#localx = np.zeros(sizex)
#localy = np.zeros(sizey)
localx = []
localy = []
strength = []
print('shape of dudx',dudx.shape)
print(sizex-1)
for i in range(1,sizex-1):
    for j in range(1,sizey-1):
        #print(i,j)
        dudx[i,j] = 0.5*(velU[i+1, j] - velU[i-1,j])/(dx[i]-dx[i-1])
        #print(i,j,dudx[i,j])
        dudy[i,j] = 0.5*(velU[i, j+1] - velU[i,j-1])/(dy[j]-dy[j-1])
        dvdx[i,j] = 0.5*(velV[i+1, j] - velV[i-1,j])/(dx[i]-dx[i-1])
        dvdy[i,j] = 0.5*(velV[i, j+1] - velV[i,j-1])/(dy[j]-dy[j-1])

#---- VORTICITY ----#
vorticity = dudy - dvdx

#A = [dudx(k) dudy(k) dudz(k); dvdx(k) dvdy(k) dvdz(k); dwdx(k) dwdy(k) -dudx(k)-dvdy(k)];

#print(dudx[10,10],dudx2[10,10])

#---- VELOCITY TENSOR GRADIENT ----#

size = np.size(velU)
for i in range(sizex):
        for j in range(sizey):
            A = [[dudx[i,j],dudy[i,j],dudz[i,j]],[dvdx[i,j],dvdy[i,j],dvdz[i,j]],[dwdx[i,j],dwdy[i,j],-dudx[i,j]-dvdy[i,j]]]
            ls = np.linalg.eigvals(A)
            ls2[i,j] = max(abs(ls.imag))
            #if ls2[i,j] > 1.0:
            #    localx = np.insert(localx,0,i)
            #    localy = np.insert(localy,0,j)
            #    strength = np.insert(strength,0,max(abs(ls.imag))) 
            #print(i,j,ls)
            #print('imag',ls2)
            #L[i,j] = 

#print('str',strength)
#print(localx)
#print(localy)
#for i in range(sizex):
#    for j in range(sizey):
#        ls2Max = np.argmax(ls2[i-1,j-1],ls2[i,j-1],ls2[i+1,j-1],ls2[i-1,j],ls2[i,j],ls2[i+1,j],ls2[i-1,j+1],ls2[i,j+1],ls2[i+1,j+1])
#        tempValue = ls2

maxima1 = ndimage.maximum_filter(ls2,size=(6,6))
maxima = find_peaks(ls2, 1.0, box_size=6)
strength = maxima['peak_value']*10
#print(maxima)
#---- PLOTTING ----#
# Make plot with vertical (default) colorbar
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
cax = ax1.imshow(velU, interpolation='nearest', cmap=plt.cm.coolwarm)
ax1.set_title('Velocity U (velocity_s)')
#cbar = fig.colorbar(cax)#, ticks=[-1, 0, 1])

cax = ax2.imshow(velV, interpolation='nearest', cmap=plt.cm.coolwarm)
ax2.set_title('Velocity V (velocity_n)')
#cbar = fig.colorbar(cax)#, ticks=[-1, 0, 1])

#cax = ax3.imshow(vorticity, interpolation='nearest', cmap=plt.cm.coolwarm)
ax3.set_title('Swirling Strength, max filter')
ax3.set_xlim(0,1152)
ax3.imshow(maxima1)
ax3.scatter(maxima['x_peak'], maxima['y_peak'],s=strength,edgecolors='r',facecolors='none')

cax = ax4.imshow(ls2, interpolation='nearest', cmap="Greys")
ax4.set_title('Swirling Strength, no filter')
#ax4.scatter(localy,localx,s=strength,edgecolors='r',facecolors='none')
#ax4.plot(maxima['x_peak'], maxima['y_peak'], ls='none',marker='o')
ax4.scatter(maxima['x_peak'], maxima['y_peak'],s=maxima['peak_value'],edgecolors='r',facecolors='none')

plt.tight_layout()
plt.show()
