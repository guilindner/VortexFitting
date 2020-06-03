#!/usr/bin/env/ python
"""class"""
import numpy as np
from netCDF4 import Dataset
import re

class VelocityField():
    """Data file

    Loads the input file with the NetCFD (.nc) format or a Tecplot format (.dat); 
    initialize the variables.
    
    :param path: file path
    :type path: str
    :param time: current time step
    :type time: int
    :param meanfilepath: in case of a mean field subtraction
    :type meanfilepath: str
    :param filetype: 'piv_netcdf', 'dns, 'dns2', 'piv_tecplot'
    :type filetype: str
    :param u: 1st component velocity field
    :type u: array of float
    :param v: 2nd component velocity field
    :type v: array of float
    :param w: 3rd component velocity field, optional
    :type w: array of float
    :param dx: spatial mesh
    :type dx: array of float
    :param dy: spatial mesh
    :type dy: array of float
    :param dz: spatial mesh, optional
    :type dz: array of float
    :param norm: for normalization of the swirling field
    :type norm: boolean
    :param normdir: False, 'x' or 'y'
    :type normdir: str
    :param step_dx: for homogeneous mesh, provides a unique step
    :type step_dx: float
    :param step_dy: for homogeneous mesh, provides a unique step 
    :type step_dy: float
    :param derivative: 
    :type derivative: array of float

    """
    def __init__(self, path="/", time=0, meanfilepath="/",filetype="/"):


        if ('{:' in path):
            self.path = path.format(time)
        else:
            self.path = path

        self.time = time
        self.meanfilepath = meanfilepath
        #filetype = 'tecplot' #change here to the desired format
        
        ## To read data
        #grp1 = Dataset(path, 'r')
        #self.u = np.array(grp1.variables['velocity_x'][time, :, :])
        #self.v = np.array(grp1.variables['velocity_y'][time, :, :])
        #self.w = np.array(grp1.variables['velocity_z'][time, :, :])

        ## To read statistics
        #grp1 = Dataset('statistics.nc', 'r')
        #self.mean = np.array(grp1.variables['um'][time, :, :])
        #grp1.close()

        ## Create a mesh
        #self.samples = self.u.shape[1]
        #self.dx = np.linspace(0, self.u.shape[1], self.u.shape[1])
        #self.dy = np.linspace(0, self.u.shape[0], self.u.shape[0])

        ## True to normalize
        #self.ifnorm = False
        #if self.ifnorm:
        #    self.normdir = False
         

        if filetype == 'piv_netcdf':
        #PIV DATA with netCDF format
            grp1 = Dataset(self.path, 'r')
            self.u = np.array(grp1.variables['velocity_n'][time, :, :])
            self.v = np.array(grp1.variables['velocity_s'][time, :, :])
            self.w = np.array(grp1.variables['velocity_z'][time, :, :])
            self.dx = np.array(grp1.variables['grid_n'])
            self.dy = np.array(grp1.variables['grid_z'])
            self.dy = self.dy - self.dy[0] #it does not start at 0
            self.u = self.u - np.mean(self.u, 1)[:, None]
            self.v = self.v - np.mean(self.v, 1)[:, None]
            self.w = self.w - np.mean(self.w, 1)[:, None]
            self.norm = True
            self.normdir = 'y'
            self.samples = self.u.shape[1]
            grp1.close()

        if filetype == 'dns':
        #DNS DATA
            grp1 = Dataset(self.path, 'r')
            self.u = np.array(grp1.variables['velocity_x'][time, :, :])
            self.v = np.array(grp1.variables['velocity_y'][time, :, :])
            self.w = np.array(grp1.variables['velocity_z'][time, :, :])
            self.samples = self.u.shape[1]
            self.dx = np.linspace(0, self.samples, self.samples)
            self.dy = np.linspace(0, self.samples, self.samples)
            self.norm = False
            self.normdir = False
            grp1.close()
            
        if filetype == 'dns2':
        # ILKAY DATA FOR DNS
            grp1 = Dataset(self.path, 'r')
            grp2 = Dataset('../data/DNS_example/vel_v_00000000.00400000.nc', 'r')
            grp3 = Dataset('../data/DNS_example/vel_w_00000000.00400000.nc', 'r')
            grp4 = Dataset('../data/DNS_example/grid_x.nc', 'r')
            grp5 = Dataset('../data/DNS_example/grid_y.nc', 'r')
            grp6 = Dataset('../data/DNS_example/grid_z.nc', 'r')
            self.u = np.array(grp1.variables['U'][60])
            self.v = np.array(grp2.variables['V'][60])
            self.w = np.array(grp3.variables['W'][60])
            self.dx = np.array(grp4.variables['gridx'][0, 0, :])
            self.dy = np.array(grp5.variables['gridy'][0, :, 0])
            self.dz = np.array(grp6.variables['gridz'][:, 0, 0])
            self.norm = False
            grp1.close()
                    
        if filetype == 'piv_tecplot':
        # DATA FORMAT FOR PIV - TECPLOT
            with open(self.path) as myfile:
                myfile.readline()
                list_variables=re.findall(r'\"(.*?)\"',myfile.readline())
		# Default: data are x, y, z, u, v, w
                index_x,index_y,index_z,index_u,index_v,index_w = 0,1,2,3,4,5 
                #for j in range(0,len(list_variables)):
                #    if list_variables[j] in ['x', 'X', 'x/c']:
                #        index_x=j
                #    if list_variables[j] in ['y', 'Y', 'y/c']:
                #        index_y=j
                #    if list_variables[j] in ['z', 'Z', 'z/c']:
                #        index_z=j                    
                #    if list_variables[j] in ['VX', 'U', 'u\'/Udeb']:
                #        index_u=j
                #    if list_variables[j] in ['VY', 'V', 'v\'/Udeb']:
                #        index_v=j
                #    if list_variables[j] in ['VZ', 'W', 'w\'/Udeb']:
                #        index_w=j                 
    
            grp1=np.loadtxt(self.path,delimiter=" ",dtype=float,skiprows=3) #skip header, default is 3 lines
            dx_tmp = np.array(grp1[:,0])
        
            for i in range(1,dx_tmp.shape[0]):
                if (dx_tmp[i]==dx_tmp[0]):
                    self.sizey=i;
                    break;
            self.sizex=np.int(dx_tmp.shape[0]/self.sizey); #domain size

            self.u  = np.array(grp1[:,index_u]).reshape(self.sizex,self.sizey)
            self.v  = np.array(grp1[:,index_v]).reshape(self.sizex,self.sizey)

            if (self.meanfilepath != '/' ):
                print("subtracting mean file")
                grp2=np.loadtxt(meanfilepath,delimiter=" ",dtype=float,skiprows=3) #mean data
                self.uMean  = np.array(grp2[:,index_u]).reshape(self.sizex,self.sizey)
                self.vMean  = np.array(grp2[:,index_v]).reshape(self.sizex,self.sizey)
                self.u = self.u - self.uMean
                self.v = self.v - self.vMean

            self.samples = self.u.shape[1]

            tmp_x  = np.array(grp1[:,index_x]).reshape(self.sizex,self.sizey)
            tmp_y  = np.array(grp1[:,index_y]).reshape(self.sizex,self.sizey)

            self.dx = np.linspace(0, np.max(tmp_x)-np.min(tmp_x), self.u.shape[1])
            self.dy = np.linspace(0, np.max(tmp_y)-np.min(tmp_y), self.u.shape[0])

            self.norm = False
            self.normdir = 'x'
     
        
        if filetype == 'test':
            sizex, sizey=100,100
            u = np.linspace(0.0,1.0,sizex)
            v = np.linspace(0.0,1.0,sizey)
            self.u, self.v = np.meshgrid(u,v)
            self.dx = np.linspace(0.0,1.0,sizex)
            self.dy = np.linspace(0.0,1.0,sizey)

       
        #COMMON TO ALL DATA
        self.step_dx=round((np.max(self.dx)-np.min(self.dx)) / (np.size(self.dx)-1) ,6)
        self.step_dy=round((np.max(self.dy)-np.min(self.dy)) / (np.size(self.dy)-1) ,6)

        self.derivative = {'dudx': np.zeros_like(self.u),
                           'dudy': np.zeros_like(self.u),
                           'dudz': np.zeros_like(self.u),
                           'dvdx': np.zeros_like(self.u),
                           'dvdy': np.zeros_like(self.u),
                           'dvdz': np.zeros_like(self.u),
                           'dwdx': np.zeros_like(self.u),
                           'dwdy': np.zeros_like(self.u),
                           'dwdz': np.zeros_like(self.u)}





