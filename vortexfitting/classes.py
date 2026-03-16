#!/usr/bin/env/ python3
"""
class VelocityField

.. tip:: new file type reader should be added here
"""

import numpy as np
import netCDF4
import h5py


class VelocityField:
    """
    Data file

    Loads the input file with the NetCFD (.nc)
    format, a Tecplot format (.dat) or a raw / dat file;
    initialize the variables.

    :param file_path: file path
    :type  file_path: str
    :param time_step: current time step
    :type  time_step: int
    :param mean_file_path: in case of a mean field subtraction
    :type  mean_file_path: str
    :param file_type: 'PIV NetCDF (.nc)', 'DNS (.nc), 'PIV Tecplot (.dat)', 'OpenFOAM (.raw)', 'HDF5 (.h5)'
    :type  file_type: str
    :param self.x_coordinate_matrix: spatial mesh
    :type  self.x_coordinate_matrix: ndarray
    :param self.y_coordinate_matrix: spatial mesh
    :type  self.y_coordinate_matrix: ndarray
    :param self.z_coordinate_matrix: spatial mesh, optional
    :type  self.z_coordinate_matrix: ndarray
    :param self.u_velocity_matrix: 1st component velocity field
    :type  self.u_velocity_matrix: ndarray
    :param self.v_velocity_matrix: 2nd component velocity field
    :type  self.v_velocity_matrix: ndarray
    :param self.w_velocity_matrix: 3rd component velocity field, optional
    :type  self.w_velocity_matrix: ndarray
    :param self.detection_field: ndarray
    :type  self.detection_field: ndarray
    :param self.normalization_flag: for normalization of the swirling field
    :type  self.normalization_flag: boolean
    :param self.normalization_direction: 'None', 'x' or 'y'
    :type  self.normalization_direction: str
    :param self.x_coordinate_step: for homogeneous mesh, provides a unique step
    :type  self.x_coordinate_step: float
    :param self.y_coordinate_step: for homogeneous mesh, provides a unique step
    :type  self.y_coordinate_step: float
    :param self.z_coordinate_step: for homogeneous mesh, provides a unique step
    :type  self.z_coordinate_step: float
    :param self.derivative: contains 'dudx', 'dudy', 'dvdx', 'dvdy'.
                       Can be extended to the 3rd dimension
    :type  self.derivative: dict

    :returns: vfield, an instance of the VelocityField class
    :rtype: class VelocityField
    """

    def __init__(self, file_path="/", time_step=0, mean_file_path="/", file_type="/"):

        if '{:' in file_path:
            self.file_path = file_path.format(time_step)
        else:
            self.file_path = file_path

        self.time_step = time_step
        self.mean_file_path = mean_file_path
        datafile_read = None
        # file_type = 'tecplot' #change here to the desired format

        # To read data
        # datafile_read = netCDF4.Dataset(path, 'r')
        # self.u_velocity_matrix = np.array(datafile_read.variables['velocity_x'][time, :, :])
        # self.v_velocity_matrix = np.array(datafile_read.variables['velocity_y'][time, :, :])
        # self.w_velocity_matrix = np.array(datafile_read.variables['velocity_z'][time, :, :])

        # To read statistics
        # datafile_read = netCDF4.Dataset('statistics.nc', 'r')
        # self.mean = np.array(datafile_read.variables['um'][time, :, :])
        # datafile_read.close()

        # Create a mesh
        # self.x_coordinate_matrix = np.linspace(0, self.u_velocity_matrix.shape[1], self.u_velocity_matrix.shape[1])
        # self.y_coordinate_matrix = np.linspace(0, self.u_velocity_matrix.shape[0], self.u_velocity_matrix.shape[0])

        # True to normalize
        # self.normalization_flag = False
        # if self.normalization_flag:
        #    self.normalization_direction = False
        print(f'Reading file: {self.file_path} file type: {file_type}')

        if file_type in ['HDF5 (.h5)','h5']:
            # print("[DEBUG] Opening HDF5 file:", self.file_path)
            try:
                with h5py.File(self.file_path, 'r') as h5_file:
                    # print("[DEBUG] Available datasets:", list(h5_file.keys()))

                    try:
                        ux = np.array(h5_file['ux'])
                        uy = np.array(h5_file['uy'])
                        self.u_velocity_matrix = ux[time_step, :, :]
                        self.v_velocity_matrix = uy[time_step, :, :]
                        _, self.y_coordinate_size, self.x_coordinate_size = ux.shape
                    except (IndexError, KeyError, ValueError) as e:
                        print("reading a mean file")
                        ux = np.array(h5_file['ux'])
                        uy = np.array(h5_file['uy'])
                        self.u_velocity_matrix = ux[:, :]
                        self.v_velocity_matrix = uy[:, :]
                        self.y_coordinate_size, self.x_coordinate_size = ux.shape
                        # print(f"[DEBUG] Loaded 'ux', shape: {h5_file['ux']}")
                        # print(f"[DEBUG] Loaded 'uy', shape: {h5_file['uy']}")
                    # self.w_velocity_matrix = uz.mean(axis=0)  # uncomment if 'uz' exists

                    # Shape

                    self.z_coordinate_size = 1
                    print("self.y_coordinate_size, self.x_coordinate_size", self.y_coordinate_size, self.x_coordinate_size)
                    # Dummy uniform coordinate mesh (structured assumption)
                    self.x_coordinate_matrix = np.linspace(0, 1, self.x_coordinate_size)
                    self.y_coordinate_matrix = np.linspace(0, 1, self.y_coordinate_size)
                    self.z_coordinate_matrix = np.zeros((self.y_coordinate_size,))

                    self.normalization_flag = False
                    self.normalization_direction = 'None'
            except Exception as e:
                print(f"\n[H5 Reading Error] Could not open file {self.file_path}: {e}")

            if self.mean_file_path and self.mean_file_path != '/':
                # print("[DEBUG] Subtracting mean HDF5 file:", self.mean_file_path)
                try:
                    with h5py.File(self.mean_file_path, 'r') as mean_h5:
                        # Load mean datasets
                        mean_ux = np.array(mean_h5['ux'])
                        mean_uy = np.array(mean_h5['uy'])

                        # Reduce to 2D if needed
                        if mean_ux.ndim > 2:
                            mean_ux = mean_ux.mean(axis=0)
                        if mean_uy.ndim > 2:
                            mean_uy = mean_uy.mean(axis=0)

                        # Subtract
                        self.u_velocity_matrix -= mean_ux
                        self.v_velocity_matrix -= mean_uy

                except Exception as e:
                    print("[ERROR] Failed to subtract mean HDF5 file:", e)

        if file_type in ['PIV NetCDF (.nc)', 'piv_netcdf']:
            # PIV DATA with netCDF format
            try:
                datafile_read = netCDF4.Dataset(self.file_path, 'r')
            except IOError:
                print("\nReading error. Maybe a wrong file type?\n")
            self.u_velocity_matrix = np.array(datafile_read.variables['velocity_n'][time_step, :, :])
            self.v_velocity_matrix = np.array(datafile_read.variables['velocity_s'][time_step, :, :])
            self.w_velocity_matrix = np.array(datafile_read.variables['velocity_z'][time_step, :, :])
            self.x_coordinate_matrix = np.array(datafile_read.variables['grid_n'])
            self.y_coordinate_matrix = np.array(datafile_read.variables['grid_z'])
            self.z_coordinate_matrix = np.zeros_like(datafile_read.variables['grid_z'])
            self.y_coordinate_matrix = self.y_coordinate_matrix - self.y_coordinate_matrix[0]  # it does not start at 0
            self.u_velocity_matrix = self.u_velocity_matrix - np.mean(self.u_velocity_matrix, 1)[:, None]
            self.v_velocity_matrix = self.v_velocity_matrix - np.mean(self.v_velocity_matrix, 1)[:, None]
            self.w_velocity_matrix = self.w_velocity_matrix - np.mean(self.w_velocity_matrix, 1)[:, None]
            self.x_coordinate_size = self.u_velocity_matrix.shape[1]
            self.y_coordinate_size = self.u_velocity_matrix.shape[0]
            self.z_coordinate_size = 1
            self.normalization_flag = True
            self.normalization_direction = 'y'
            datafile_read.close()

        if file_type in ['DNS (.nc)', 'dns']:
            # DNS DATA, netCDF format
            try:
                datafile_read = netCDF4.Dataset(self.file_path, 'r')
            except IOError:
                print("\nReading error. Maybe a wrong file type?\n")

            self.u_velocity_matrix = np.array(datafile_read.variables['velocity_x'][time_step, :, :])
            self.v_velocity_matrix = np.array(datafile_read.variables['velocity_y'][time_step, :, :])
            self.w_velocity_matrix = np.array(datafile_read.variables['velocity_z'][time_step, :, :])
            self.x_coordinate_matrix = np.linspace(0, self.u_velocity_matrix.shape[1], self.u_velocity_matrix.shape[1])
            self.y_coordinate_matrix = np.linspace(0, self.u_velocity_matrix.shape[0], self.u_velocity_matrix.shape[0])
            self.z_coordinate_matrix = np.linspace(0, self.u_velocity_matrix.shape[0], self.u_velocity_matrix.shape[0])
            self.x_coordinate_size = self.u_velocity_matrix.shape[1]
            self.y_coordinate_size = self.u_velocity_matrix.shape[0]
            self.z_coordinate_size = 1
            self.normalization_flag = False
            self.normalization_direction = 'None'
            datafile_read.close()
            if self.mean_file_path != '/':
                print("subtracting mean file")  # load and subtract mean data
                datafile_mean_read = netCDF4.Dataset(self.mean_file_path, 'r')
                u_velocity_matrix_mean = np.array(datafile_mean_read.variables['velocity_x'][:, :])
                v_velocity_matrix_mean = np.array(datafile_mean_read.variables['velocity_y'][:, :])
                self.u_velocity_matrix = self.u_velocity_matrix - u_velocity_matrix_mean
                self.v_velocity_matrix = self.v_velocity_matrix - v_velocity_matrix_mean
                if self.w_velocity_matrix is not None:
                    w_velocity_matrix_mean = np.array(datafile_read.variables['velocity_z'][:, :])
                    self.w_velocity_matrix = self.w_velocity_matrix - w_velocity_matrix_mean

        # if file_type == 'dns2':
        #    DATA FOR DNS, netCDF format
        #    datafile_read = netCDF4.Dataset(self.file_path, 'r')
        #    grp2 = netCDF4.Dataset('../data/DNS_example/vel_v_00000000.00400000.nc', 'r')
        #    grp3 = netCDF4.Dataset('../data/DNS_example/vel_w_00000000.00400000.nc', 'r')
        #    grp4 = netCDF4.Dataset('../data/DNS_example/grid_x.nc', 'r')
        #    grp5 = netCDF4.Dataset('../data/DNS_example/grid_y.nc', 'r')
        #    grp6 = netCDF4.Dataset('../data/DNS_example/grid_z.nc', 'r')
        #    self.u_velocity_matrix = np.array(datafile_read.variables['U'][60])
        #    self.v_velocity_matrix = np.array(grp2.variables['V'][60])
        #    self.w_velocity_matrix = np.array(grp3.variables['W'][60])
        #    self.x_coordinate_matrix = np.array(grp4.variables['gridx'][0, 0, :])
        #    self.y_coordinate_matrix = np.array(grp5.variables['gridy'][0, :, 0])
        #    self.z_coordinate_matrix = np.array(grp6.variables['gridz'][:, 0, 0])
        #    self.normalization_flag = False
        #    datafile_read.close()

        if file_type in ['PIV Tecplot (.dat)','piv_tecplot']:
            # DATA FORMAT FOR PIV - TECPLOT
            # if you want to read the variables list and look automatically to the indexes
            # (please uncomment also the "import re", as you need regular expressions
            # with open(self.file_path) as myfile:
            # myfile.readline()
            # list_variables=re.findall(r'\"(.*?)\"',myfile.readline())
            # Default: data are x, y, z, u, v, w
            # index_x,index_y,index_z,index_u,index_v,index_w = 0,1,2,3,4,5
            # for j in range(0,len(list_variables)):
            #    if list_variables[j] in ['x', 'X']:
            #        index_x=j
            #    if list_variables[j] in ['y', 'Y']:
            #        index_y=j
            #    if list_variables[j] in ['z', 'Z']:
            #        index_z=j
            #    if list_variables[j] in ['VX', 'U']:
            #        index_u=j
            #    if list_variables[j] in ['VY', 'V']:
            #        index_v=j
            #    if list_variables[j] in ['VZ', 'W']:
            #        index_w=j
            try:
                datafile_read = np.loadtxt(
                    self.file_path, delimiter=" ", dtype=float, skiprows=3
                )  # skip header, default is 3 lines
            except IOError:
                print("\nReading error. Maybe a wrong file type?\n")

            index_x, index_y, index_z, index_u, index_v, index_w = 0, 1, 2, 3, 4, 5
            dx_tmp = np.array(datafile_read[:, index_x])

            for i in range(1, dx_tmp.shape[0]):
                if dx_tmp[i] == dx_tmp[0]:
                    self.y_coordinate_size = i
                    break
            self.x_coordinate_size = int(dx_tmp.shape[0] / self.y_coordinate_size)  # domain size
            self.z_coordinate_size = 1

            self.u_velocity_matrix = np.array(datafile_read[:, index_u]).reshape(
                self.x_coordinate_size, self.y_coordinate_size
            )
            self.v_velocity_matrix = np.array(datafile_read[:, index_v]).reshape(
                self.x_coordinate_size, self.y_coordinate_size
            )
            try:
                self.w_velocity_matrix = np.array(datafile_read[:, index_w]).reshape(
                    self.x_coordinate_size, self.y_coordinate_size
                )
            except IndexError:
                print('No w velocity matrix')

            if self.mean_file_path != '/':
                print("subtracting mean file")
                # load and subtract mean data
                datafile_mean_read = np.loadtxt(mean_file_path, delimiter=" ", dtype=float, skiprows=3)
                u_velocity_matrix_mean = np.array(datafile_mean_read[:, index_u]).reshape(
                    self.x_coordinate_size, self.y_coordinate_size
                )
                v_velocity_matrix_mean = np.array(datafile_mean_read[:, index_v]).reshape(
                    self.x_coordinate_size, self.y_coordinate_size
                )
                self.u_velocity_matrix = self.u_velocity_matrix - u_velocity_matrix_mean
                self.v_velocity_matrix = self.v_velocity_matrix - v_velocity_matrix_mean
                if self.w_velocity_matrix is not None:
                    w_velocity_matrix_mean = np.array(datafile_mean_read[:, index_w]).reshape(
                        self.x_coordinate_size, self.y_coordinate_size
                    )
                    self.w_velocity_matrix = self.w_velocity_matrix - w_velocity_matrix_mean

            tmp_x = np.array(datafile_read[:, index_x]).reshape(self.x_coordinate_size, self.y_coordinate_size)
            tmp_y = np.array(datafile_read[:, index_y]).reshape(self.x_coordinate_size, self.y_coordinate_size)
            self.x_coordinate_matrix = np.linspace(0, np.max(tmp_x) - np.min(tmp_x), self.u_velocity_matrix.shape[1])
            self.y_coordinate_matrix = np.linspace(0, np.max(tmp_y) - np.min(tmp_y), self.u_velocity_matrix.shape[0])
            try:
                tmp_z = np.array(datafile_read[:, index_z]).reshape(self.x_coordinate_size, self.y_coordinate_size)
                self.z_coordinate_matrix = tmp_z[0, 0]
            except IndexError:
                print('No z component')

            self.normalization_flag = False
            self.normalization_direction = 'None'

        if file_type in ['OpenFOAM (.raw)', 'openfoam']:
            try:
                datafile_read = np.loadtxt(
                    self.file_path, delimiter=" ", dtype=float, skiprows=2
                )  # skip header, default is 2 lines
            except IOError:
                print("\nReading error. Maybe a wrong file type?\n")

            index_x, index_y, index_z, index_u, index_v, index_w = 0, 1, 2, 3, 4, 5
            dx_tmp = np.array(datafile_read[:, index_x])
            for i in range(1, dx_tmp.shape[0]):
                if dx_tmp[i] == dx_tmp[0]:
                    self.y_coordinate_size = i
                    break
            self.x_coordinate_size = int(dx_tmp.shape[0] / self.y_coordinate_size)  # domain size
            self.z_coordinate_size = 1

            self.u_velocity_matrix = np.array(datafile_read[:, index_u]).reshape(
                self.x_coordinate_size, self.y_coordinate_size
            )
            self.v_velocity_matrix = np.array(datafile_read[:, index_v]).reshape(
                self.x_coordinate_size, self.y_coordinate_size
            )

            try:
                self.w_velocity_matrix = np.array(datafile_read[:, index_w]).reshape(
                    self.x_coordinate_size, self.y_coordinate_size
                )
            except IndexError:
                print('No w velocity matrix')
            tmp_x = np.array(datafile_read[:, index_x]).reshape(self.x_coordinate_size, self.y_coordinate_size)
            tmp_y = np.array(datafile_read[:, index_y]).reshape(self.x_coordinate_size, self.y_coordinate_size)
            self.x_coordinate_matrix = np.linspace(0, np.max(tmp_x) - np.min(tmp_x), self.u_velocity_matrix.shape[1])
            self.y_coordinate_matrix = np.linspace(0, np.max(tmp_y) - np.min(tmp_y), self.u_velocity_matrix.shape[0])
            try:
                tmp_z = np.array(datafile_read[:, index_z]).reshape(self.x_coordinate_size, self.y_coordinate_size)
                self.z_coordinate_matrix = tmp_z[0, 0]
            except IndexError:
                print('No z component')

            self.normalization_flag = False
            self.normalization_direction = 'None'

        if file_type == 'test':
            self.x_coordinate_size, self.y_coordinate_size, self.z_coordinate_size = 11, 5, 1
            u = np.linspace(0.0, 1.0, self.x_coordinate_size)
            v = np.linspace(0.0, 1.0, self.y_coordinate_size)
            self.u_velocity_matrix, self.v_velocity_matrix = np.meshgrid(u, v)
            self.x_coordinate_matrix = np.linspace(0.0, 1.0, self.x_coordinate_size)
            self.y_coordinate_matrix = np.linspace(0.0, 1.0, self.y_coordinate_size)

        # COMMON TO ALL DATA
        self.x_coordinate_step = round(
            (np.max(self.x_coordinate_matrix).astype(float) - np.min(self.x_coordinate_matrix).astype(float))
            / (np.size(self.x_coordinate_matrix) - 1),
            6,
        )
        self.y_coordinate_step = round(
            (np.max(self.y_coordinate_matrix).astype(float) - np.min(self.y_coordinate_matrix).astype(float))
            / (np.size(self.y_coordinate_matrix) - 1),
            6,
        )
        self.z_coordinate_step = 0.0

        self.derivative = {
            'dudx': np.zeros_like(self.u_velocity_matrix),
            'dudy': np.zeros_like(self.u_velocity_matrix),
            'dudz': np.zeros_like(self.u_velocity_matrix),
            'dvdx': np.zeros_like(self.u_velocity_matrix),
            'dvdy': np.zeros_like(self.u_velocity_matrix),
            'dvdz': np.zeros_like(self.u_velocity_matrix),
            'dwdx': np.zeros_like(self.u_velocity_matrix),
            'dwdy': np.zeros_like(self.u_velocity_matrix),
            'dwdz': np.zeros_like(self.u_velocity_matrix),
        }
