from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

import fitting

grp1 = Dataset('../data/generatedField.nc', 'w', format='NETCDF4')
grp1.description = 'Sample field with an oseen vortex'

ndim = 256 # spacing

# dimensions
grp1.createDimension('resolution_x', ndim)
grp1.createDimension('resolution_y', ndim)
grp1.createDimension('resolution_z', 1)

# variables
velocity_x = grp1.createVariable('velocity_x', 'f4', ('resolution_z','resolution_y','resolution_x'))
velocity_y = grp1.createVariable('velocity_y', 'f4', ('resolution_z','resolution_y','resolution_x'))
velocity_z = grp1.createVariable('velocity_z', 'f4', ('resolution_z','resolution_y','resolution_x'))

# data
velocity_x[:] = np.random.random((1,ndim,ndim))/10
velocity_y[:] = np.random.random((1,ndim,ndim))/10
velocity_z[:] = np.random.random((1,ndim,ndim))/10

# grid
x = np.linspace(0,ndim,ndim)
y = np.linspace(0,ndim,ndim)

#dist = 40
#X = np.linspace(-1,1,dist)
#Y = np.linspace(-1,1,dist)
xx, yy = np.meshgrid(x,y)
coreR = 1.0
gamma = 50
fxCenter = 64
fyCenter = 192
u_conv = 0.0
v_conv = 0.0
Uw, Vw = fitting.velocity_model(coreR, gamma, fxCenter, fyCenter, u_conv, v_conv, xx, yy)
Uw = Uw + u_conv
Vw = Vw + v_conv
#xCenter = 200 #where to move the vortex
#yCenter = 100
#print(Uw)
#print(Vw)
#for i in range(dist):
#    for j in range(dist):
#        velocity_x[0,i+xCenter,j+yCenter] = Uw[i,j]
#        velocity_y[0,i+xCenter,j+yCenter] = Vw[i,j]
#x = np.linspace(0,ndim,ndim)
#y = np.linspace(0,ndim,ndim)
#xx, yy = np.meshgrid(x,y)
#velx = velocity_x[0]
#vely = velocity_y[0]
velocity_x[0,:,:] += Uw[:,:]
velocity_y[0,:,:] += Vw[:,:]
s = 4
#velx = np.einsum('ij->ji',velx)
#vely = np.einsum('ij->ji',vely)
plt.quiver(xx[::s,::s],yy[::s,::s],velocity_x[0,::s,::s],velocity_y[0,::s,::s])
#plt.quiver(X,Y,Uw,Vw)
plt.show()
grp1.close()
