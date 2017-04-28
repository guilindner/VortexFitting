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
            self.u = np.array(grp1.variables['velocity_s'][self.time])
            self.v = np.array(grp1.variables['velocity_n'][self.time])
            self.u = np.einsum('ij->ji',self.u)
            self.v = np.einsum('ij->ji',self.v)
            self.dx = np.array(grp1.variables['grid_n'])
            self.dy = np.array(grp1.variables['grid_z'])
            self.sizex = self.u.shape[0]
            self.sizey = self.u.shape[1]
            
        elif 'U' in grp1.variables.keys():
            grp2 = Dataset('../data/DNS_example/vel_v_00000000.00400000.nc','r')
            grp3 = Dataset('../data/DNS_example/vel_w_00000000.00400000.nc','r')
            grp4 = Dataset('../data/DNS_example/grid_x.nc','r')
            grp5 = Dataset('../data/DNS_example/grid_y.nc','r')
            grp6 = Dataset('../data/DNS_example/grid_z.nc','r')
            self.u = np.array(grp1.variables['U'][80])
            self.v = np.array(grp2.variables['V'][80])
            self.w = np.array(grp3.variables['W'][80])
            domain_x = 1.
            domain_y = 1.
            domain_z = 1.
            self.dx = np.linspace(0,domain_x,grp1.dimensions['resolution_x'].size)
            self.dy = np.linspace(0,domain_y,grp1.dimensions['resolution_y'].size)
            self.dz = np.linspace(0,domain_z,grp1.dimensions['resolution_z'].size)
            self.sizex = self.u.shape[0]
            self.sizey = self.u.shape[1]
            self.sizez = self.u.shape[1]#fix later for different size
        else:
            print('Netcdf file format not recognized')		
        
        self.derivative = {'dudx': np.zeros((self.sizex,self.sizey)),
                           'dudy': np.zeros((self.sizex,self.sizey)),
                           'dudz': np.zeros((self.sizex,self.sizey)),
                           'dvdx': np.zeros((self.sizex,self.sizey)),
                           'dvdy': np.zeros((self.sizex,self.sizey)),
                           'dvdz': np.zeros((self.sizex,self.sizey)),
                           'dwdx': np.zeros((self.sizex,self.sizey)),
                           'dwdy': np.zeros((self.sizex,self.sizey)),
                           'dwdz': np.zeros((self.sizex,self.sizey))}
