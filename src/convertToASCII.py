import numpy as np
from netCDF4 import Dataset

grp1 = Dataset("test_data.nc",'r') 

samples = grp1.variables['velocity_x'].shape[0]
dx = np.linspace(0,10,self.samples)
dy = np.linspace(0,10,self.samples)
u = np.array(grp1.variables['velocity_x'][:,:,:])
v = np.array(grp1.variables['velocity_y'][:,:,:])
w = np.array(grp1.variables['velocity_z'][:,:,:])
u = u[0]
v = v[0]
w = w[0]
u = np.einsum('ij->ji',self.u)
v = np.einsum('ij->ji',self.v)

file = open('ascii','w')
