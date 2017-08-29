#!/usr/bin/env/ python
"""class"""
import numpy as np
from netCDF4 import Dataset

class VelocityField():
    """NetCDF file 
    
    Loads the input file with the NetCFD (.nc) format and
    initialize the variables.
    
    """
    def __init__(self,path="/",time=0):
        self.path = path
        self.time = time
        grp1 = Dataset(path,'r') 
        
        self.u = np.array(grp1.variables['velocity_x'][time,:,:])
        self.v = np.array(grp1.variables['velocity_y'][time,:,:])
        self.w = np.array(grp1.variables['velocity_z'][time,:,:])
        self.samples = self.u.shape[1]
        self.dx = np.linspace(0,self.samples,self.samples)
        self.dy = np.linspace(0,self.samples,self.samples)
        self.norm = False
        self.normdir = 1

        ##EXAMPLES
        #PIV DATA   
        #self.u = np.array(grp1.variables['velocity_n'][time,:,:])
        #self.v = np.array(grp1.variables['velocity_s'][time,:,:])
        #self.w = np.array(grp1.variables['velocity_z'][time,:,:])
        #self.dx = np.array(grp1.variables['grid_n'])
        #self.dy = np.array(grp1.variables['grid_z'])
        #self.dy = self.dy - self.dy[0] #it does not start at 0
        #self.u = self.u - np.mean(self.u,1)[:,None]
        #self.v = self.v - np.mean(self.v,1)[:,None]
        #self.w = self.w - np.mean(self.w,1)[:,None]
        #self.norm = True
        #self.normdir = 0
        #self.samples = self.u.shape[1]
        
        #DNS DATA  
        #self.u = np.array(grp1.variables['velocity_x'][time,:,:])
        #self.v = np.array(grp1.variables['velocity_y'][time,:,:])
        #self.w = np.array(grp1.variables['velocity_z'][time,:,:])
        #self.samples = self.u.shape[1]
        #self.dx = np.linspace(0,self.samples,self.samples)
        #self.dy = np.linspace(0,self.samples,self.samples)
        #self.norm = False
        #self.normdir = 1
        
        # ILKAY DATA        
        #grp2 = Dataset('../data/DNS_example/vel_v_00000000.00400000.nc','r')
        #grp3 = Dataset('../data/DNS_example/vel_w_00000000.00400000.nc','r')
        #grp4 = Dataset('../data/DNS_example/grid_x.nc','r')
        #grp5 = Dataset('../data/DNS_example/grid_y.nc','r')
        #grp6 = Dataset('../data/DNS_example/grid_z.nc','r')
        #self.u = np.array(grp1.variables['U'][60])
        #self.v = np.array(grp2.variables['V'][60])
        #self.w = np.array(grp3.variables['W'][60])
        #self.dx = np.array(grp4.variables['gridx'][0,0,:])
        #self.dy = np.array(grp5.variables['gridy'][0,:,0])
        #self.dz = np.array(grp6.variables['gridz'][:,0,0])
        #self.norm = False
        
        
        self.derivative = {'dudx': np.zeros_like(self.u),
                           'dudy': np.zeros_like(self.u),
                           'dudz': np.zeros_like(self.u),
                           'dvdx': np.zeros_like(self.u),
                           'dvdy': np.zeros_like(self.u),
                           'dvdz': np.zeros_like(self.u),
                           'dwdx': np.zeros_like(self.u),
                           'dwdy': np.zeros_like(self.u),
                           'dwdz': np.zeros_like(self.u)}
