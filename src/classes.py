#!/usr/bin/env/ python
"""class"""
import numpy as np
from netCDF4 import Dataset

class VelocityField():
    """NetCDF file 
    
    Loads the input file with the NetCFD (.nc) format and
    initialize the variables.
    
    TODO: Detect PIV or DNS data.
    
    """
    def __init__(self,path="/",time=0):
        self.path = path
        self.time = time
        grp1 = Dataset(path,'r') 

        if 'velocity_s' in grp1.variables.keys():
            
            self.samples = grp1.variables['velocity_s'].shape[0]
            self.dx = np.array(grp1.variables['grid_z'])
            self.dy = np.array(grp1.variables['grid_n'])
            self.u = np.zeros((self.dx.size,self.dy.size))
            self.v = np.zeros((self.dx.size,self.dy.size))
            self.z = np.zeros((self.dx.size,self.dy.size))
            u = np.array(grp1.variables['velocity_s'][:,:,:])
            v = np.array(grp1.variables['velocity_n'][:,:,:])
            w = np.array(grp1.variables['velocity_z'][:,:,:])
            u = u - np.mean(u,axis=(0,1))[None,None,:]
            v = v - np.mean(v,axis=(0,1))[None,None,:]
            w = w - np.mean(w,axis=(0,1))[None,None,:]
            self.u = u[15]
            self.v = v[15]
            self.w = w[15]

        elif 'velocity_x' in grp1.variables.keys():
            self.samples = grp1.variables['velocity_x'].shape[0]
            self.dx = np.linspace(0,10,self.samples)
            self.dy = np.linspace(0,10,self.samples)
            u = np.array(grp1.variables['velocity_x'][:,:,:])
            v = np.array(grp1.variables['velocity_y'][:,:,:])
            w = np.array(grp1.variables['velocity_z'][:,:,:])
            #u = u - np.mean(u,axis=(0,1))[None,None,:]
            #v = v - np.mean(v,axis=(0,1))[None,None,:]
            #w = w - np.mean(w,axis=(0,1))[None,None,:]
            self.u = u[0]
            self.v = v[0]
            self.w = w[0]
            self.u = np.einsum('ij->ji',self.u)
            self.v = np.einsum('ij->ji',self.v)
                
        elif 'U' in grp1.variables.keys():
            grp2 = Dataset('../data/DNS_example/vel_v_00000000.00400000.nc','r')
            grp3 = Dataset('../data/DNS_example/vel_w_00000000.00400000.nc','r')
            grp4 = Dataset('../data/DNS_example/grid_x.nc','r')
            grp5 = Dataset('../data/DNS_example/grid_y.nc','r')
            grp6 = Dataset('../data/DNS_example/grid_z.nc','r')
            self.u = np.array(grp1.variables['U'][60])
            self.v = np.array(grp2.variables['V'][60])
            self.w = np.array(grp3.variables['W'][60])
            self.dx = np.array(grp4.variables['gridx'][0,0,:])
            self.dy = np.array(grp5.variables['gridy'][0,:,0])
            self.dz = np.array(grp6.variables['gridz'][:,0,0])

        else:
            print('Netcdf file format not recognized')		
        
        self.derivative = {'dudx': np.zeros_like(self.u),
                           'dudy': np.zeros_like(self.u),
                           'dudz': np.zeros_like(self.u),
                           'dvdx': np.zeros_like(self.u),
                           'dvdy': np.zeros_like(self.u),
                           'dvdz': np.zeros_like(self.u),
                           'dwdx': np.zeros_like(self.u),
                           'dwdy': np.zeros_like(self.u),
                           'dwdz': np.zeros_like(self.u)}
