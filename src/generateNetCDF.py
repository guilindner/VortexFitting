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
velocity_x[:] = 0#np.random.random((1,ndim,ndim))
velocity_y[:] = 0#np.random.random((1,ndim,ndim))
velocity_z[:] = np.random.random((1,ndim,ndim))

# insert oseen vortex
dist = 12
X = np.linspace(-1,1,dist)
Y = np.linspace(-1,1,dist)
X, Y = np.meshgrid(X,Y)
coreR = 1
gamma = -30
fxCenter = 0.0
fyCenter = 0.0
u_conv = 0.0
v_conv = 0.0
Uw, Vw = fitting.velocity_model(coreR, gamma, fxCenter, fyCenter, u_conv, v_conv, X, Y)
Uw = Uw + u_conv
Vw = Vw + v_conv
xCenter = 200 #where to move the vortex
yCenter = 100
#print(Uw)
#print(Vw)
for i in range(dist):
    for j in range(dist):
        velocity_x[0,i+xCenter,j+yCenter] = Uw[i,j]
        velocity_y[0,i+xCenter,j+yCenter] = Vw[i,j]
x = np.linspace(0,ndim,ndim)
y = np.linspace(0,ndim,ndim)
xx, yy = np.meshgrid(x,y)
velx = velocity_x[0]
vely = velocity_y[0]
s = 4
velx = np.einsum('ij->ji',velx)
vely = np.einsum('ij->ji',vely)
plt.quiver(xx[::s,::s],yy[::s,::s],vely[::s,::s],velx[::s,::s])
#plt.quiver(Y,X,Vw,Uw)
plt.show()
grp1.close()
