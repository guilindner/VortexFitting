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
velocity_x[:] = np.random.random((1,ndim,ndim))/1
velocity_y[:] = np.random.random((1,ndim,ndim))/1
velocity_z[:] = np.random.random((1,ndim,ndim))/1

# grid
x = np.linspace(0,ndim,ndim)
y = np.linspace(0,ndim,ndim)

#dist = 40
#x_index = np.linspace(-1,1,dist)
#y_index = np.linspace(-1,1,dist)
xx, yy = np.meshgrid(x,y)
coreR = 5.0
gamma = 30
x_real = 64
y_real = 192
u_conv = 0.5
v_conv = 0.0
u_data, v_data = fitting.velocity_model(coreR, gamma, x_real, y_real, u_conv, v_conv, xx, yy)
u_data = u_data + u_conv
v_data = v_data + v_conv
#x_center_index = 200 #where to move the vortex
#y_center_index = 100
#print(u_data)
#print(v_data)
#for i in range(dist):
#    for j in range(dist):
#        velocity_x[0,i+x_center_index,j+y_center_index] = u_data[i,j]
#        velocity_y[0,i+x_center_index,j+y_center_index] = v_data[i,j]
#x = np.linspace(0,ndim,ndim)
#y = np.linspace(0,ndim,ndim)
#xx, yy = np.meshgrid(x,y)
#velx = velocity_x[0]
#vely = velocity_y[0]
velocity_x[0,:,:] += u_data[:,:]
velocity_y[0,:,:] += v_data[:,:]
s = 4
#velx = np.einsum('ij->ji', velx)
#vely = np.einsum('ij->ji', vely)
plt.quiver(xx[::s,::s],yy[::s,::s],velocity_x[0,::s,::s],velocity_y[0,::s,::s])
#plt.quiver(x_index, y_index, u_data, v_data)
plt.show()
grp1.close()
