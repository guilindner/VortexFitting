from netCDF4 import Dataset
import numpy as np

grp1 = Dataset('../data/anand_data.nc', 'w', format='NETCDF4')
grp1.description = 'Experiments conducted at Rouen by Anand'

ndimx = 159 # spacing
ndimy = 134 # spacing

# dimensions
grp1.createDimension('resolution_x', ndimx)
grp1.createDimension('resolution_y', ndimy)
grp1.createDimension('resolution_z', 1)

# variables
velocity_x = grp1.createVariable('velocity_x', 'f4', ('resolution_z',
                                                      'resolution_y',
                                                      'resolution_x'))
velocity_y = grp1.createVariable('velocity_y', 'f4', ('resolution_z',
                                                      'resolution_y',
                                                      'resolution_x'))
velocity_z = grp1.createVariable('velocity_z', 'f4', ('resolution_z',
                                                      'resolution_y',
                                                      'resolution_x'))
grid_x = grp1.createVariable('grid_x', 'f4', 'resolution_x')
grid_y = grp1.createVariable('grid_y', 'f4', 'resolution_y')

# data
#velocity_x[:] = np.random.random((1,ndimy,ndimx))/1
#velocity_y[:] = np.random.random((1,ndimy,ndimx))/1
#velocity_z[:] = np.random.random((1,ndimy,ndimx))/1

# grid
x = np.linspace(0, ndimy, ndimx)
y = np.linspace(0, ndimy, ndimx)

infile = open('../data/guil_anand_data.dat', 'r')
lines = infile.readlines()
for j in range(ndimy):
    for i in range(ndimx):
        velocity_x[0, j, i] = lines[j*ndimx+i].split()[2]
        velocity_y[0, j, i] = lines[j*ndimx+i].split()[3]
        if j == 0:
            grid_x[i] = lines[i].split()[0]
        if i == 0:
            grid_y[j] = lines[j*ndimx].split()[1]

grp1.close()
