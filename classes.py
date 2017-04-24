#!/usr/bin/env/ python
"""class"""
import numpy as np
from netCDF4 import Dataset

class VelocityField():
    """NetCDF file """
    def __init__(self,path="/",time=0):
        self.path = path
        self.time = time
        grp1 = Dataset(path,'r') 
        self.u = np.array(grp1.variables['velocity_s'][self.time])
        self.v = np.array(grp1.variables['velocity_n'][self.time])
        self.u = np.einsum('ij->ji',self.u)
        self.v = np.einsum('ij->ji',self.v)
        self.dx = np.array(grp1.variables['grid_n'])
        self.dy = np.array(grp1.variables['grid_z'])
        self.sizex = self.u.shape[0]
        self.sizey = self.u.shape[1]
        
        self.derivative = {'dudx': np.zeros((self.sizex,self.sizey)),
                           'dudy': np.zeros((self.sizex,self.sizey)),
                           'dudz': np.zeros((self.sizex,self.sizey)),
                           'dvdx': np.zeros((self.sizex,self.sizey)),
                           'dvdy': np.zeros((self.sizex,self.sizey)),
                           'dvdz': np.zeros((self.sizex,self.sizey)),
                           'dwdx': np.zeros((self.sizex,self.sizey)),
                           'dwdy': np.zeros((self.sizex,self.sizey)),
                           'dwdz': np.zeros((self.sizex,self.sizey))}
