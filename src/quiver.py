import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
path = '../data/test_data.nc'
grp1 = Dataset(path,'r')
dx = np.array(grp1.variables['grid_z'])
dy = np.array(grp1.variables['grid_n'])
u = np.array(grp1.variables['velocity_s'][:,:,:])
v = np.array(grp1.variables['velocity_n'][:,:,:])
w = np.array(grp1.variables['velocity_z'][:,:,:])

u = u - np.mean(u,axis=(0,1))[None,None,:]
v = v - np.mean(v,axis=(0,1))[None,None,:]
w = w - np.mean(w,axis=(0,1))[None,None,:]

xCenter = 125
yCenter = 125
dist = 125

X, Y = np.meshgrid(dx[xCenter-dist:xCenter+dist],
                       dy[yCenter-dist:yCenter+dist])
U = u[15,xCenter-dist:xCenter+dist,yCenter-dist:yCenter+dist]
V = v[15,xCenter-dist:xCenter+dist,yCenter-dist:yCenter+dist]
plt.figure()
plt.title('Arrows scale with plot width, not view')
s = 8

Q = plt.quiver(X[::s,::s], Y[::s,::s], U[::s,::s], V[::s,::s])#, units='xy', scale_units='xy', scale=1)
 
plt.show()
