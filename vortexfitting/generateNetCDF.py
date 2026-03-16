import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
import netCDF4
import fitting
matplotlib.use('WXAgg')


def generateNetCDF(outfile: str, ndim: int, core_radius: float, gamma: float, x_center: float, y_center: float, 
                   u_advection: float, v_advection: float, uz0: float, theoretical_model: str) -> None:
    """
    Produce a NetCDF file with a theoretical model of a vortex

    Args:
    :param outfile: path for the out file
    :param ndim: grid size
    :param core_radius: vortex radius
    :param gamma: vortex circulation
    :param x_center: x coordinate of the vortex center
    :param y_center: y coordinate of the vortex center
    :param u_advection: u component for advection velocity
    :param v_advection: v component for advection velocity
    :param uz0: vertical velocity
    :param theoretical_model: chosen model, could be Rankine, Lamb-Oseen or Batchelor
    :type outfile: str
    :type ndim: int
    :type core_radius: float
    :type gamma: float
    :type x_center: float
    :type y_center: float
    :type u_advection: float
    :type v_advection: float
    :type uz0: float
    :type theoretical_model: str
               
    :returns: file
    :rtype: NetCDF file
    """
    print(f"Producing {outfile} with a {ndim}x{ndim} grid")
    
    # Création du fichier NetCDF
    try:
        datafile_write = netCDF4.Dataset(outfile, 'w', format='NETCDF4')
    except IOError as e:
        raise IOError(f"Error when writing file : {e}")
    
    datafile_write.description = 'Sample field with an Oseen vortex'
    
    # Création des dimensions
    datafile_write.createDimension('resolution_x', ndim)
    datafile_write.createDimension('resolution_y', ndim)
    datafile_write.createDimension('resolution_z', 1)
    
    # Création des variables
    velocity_x = datafile_write.createVariable('velocity_x', 'f4', ('resolution_z', 'resolution_y', 'resolution_x'))
    velocity_y = datafile_write.createVariable('velocity_y', 'f4', ('resolution_z', 'resolution_y', 'resolution_x'))
    velocity_z = datafile_write.createVariable('velocity_z', 'f4', ('resolution_z', 'resolution_y', 'resolution_x'))
    
    # Données aléatoires initiales
    velocity_x[:] = np.random.random((1, ndim, ndim)) / 10
    velocity_y[:] = np.random.random((1, ndim, ndim)) / 10
    velocity_z[:] = np.random.random((1, ndim, ndim)) / 10
    
    # Grille
    x_grid = np.linspace(0, ndim, ndim)
    y_grid = np.linspace(0, ndim, ndim)
    x_matrix, y_matrix = np.meshgrid(x_grid, y_grid)
        
    if theoretical_model == 'rankine':
        u_data, v_data = fitting.velocity_model(core_radius, gamma, x_center, y_center, u_advection, v_advection, 
                                                x_matrix, y_matrix, None, theoretical_model)
        w_data = np.zeros_like(u_data)
    elif theoretical_model == 'batchelor':
        u_data, v_data, w_data = fitting.velocity_model(core_radius, gamma, x_center, y_center, 
                                                        u_advection, v_advection, x_matrix, y_matrix, 
                                                        uz0, theoretical_model)        
    else:
        u_data, v_data = fitting.velocity_model(core_radius, gamma, x_center, y_center, u_advection, v_advection, 
                                                x_matrix, y_matrix, None, theoretical_model)            
        w_data = np.zeros_like(u_data)        
    u_data = u_data + u_advection
    v_data = v_data + v_advection
    # Ajout des données calculées
    velocity_x[0, :, :] += u_data[:, :]
    velocity_y[0, :, :] += v_data[:, :]
    velocity_z[0, :, :] += w_data[:, :]
    # Visualisation
    s = 4  # facteur d'échantillonnage pour quiver plot
    plt.contourf(x_matrix, y_matrix, velocity_z[0, :, :], cmap='jet')
    plt.quiver(x_matrix[::s, ::s], y_matrix[::s, ::s], velocity_x[0, ::s, ::s], velocity_y[0, ::s, ::s])
    plt.savefig(outfile.replace('.nc', '.png'))
    plt.close('all')
    
    # Fermeture du fichier
    datafile_write.close()
    print(f"Fichier {outfile} généré avec succès.")
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate a NetCDF file with a theoretical vortex model',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-o', '--outfile', dest='outfile', type=str,
                        help='Output NetCDF file', metavar='FILE', default='../data/test_netCDF.nc')
    parser.add_argument('-n', '--ndim', dest='ndim', type=int, default=256,
                        help='Grid size (default: 256)')
    parser.add_argument('-r', '--radius', dest='core_radius', type=float, default=5.0,
                        help='Vortex core radius (default: 5.0)')
    parser.add_argument('-g', '--gamma', dest='gamma', type=float, default=30.0,
                        help='Vortex circulation (default: 30.0)')
    parser.add_argument('-xc', '--xcenter', dest='x_center', type=float, default=64.0,
                        help='X-coordinate of vortex center (default: 64.0)')
    parser.add_argument('-yc', '--ycenter', dest='y_center', type=float, default=192.0,
                        help='Y-coordinate of vortex center (default: 192.0)')
    parser.add_argument('-ua', '--uadvection', dest='u_advection', type=float, default=0.0,
                        help='U component of advection velocity (default: 0.0)')
    parser.add_argument('-va', '--vadvection', dest='v_advection', type=float, default=0.0,
                        help='V component of advection velocity (default: 0.0)')
    parser.add_argument('-uz', '--uz0', dest='uz0', type=float, default=0.5,
                        help='Vertical velocity component (default: 0.5)')
    parser.add_argument('-m', '--model', dest='theoretical_model', type=str, 
                        choices=['rankine', 'lamb-oseen', 'batchelor'], 
                        help='Theoretical vortex model (rankine, lamb-oseen, batchelor)', 
                        default='lamb-oseen')

    args = parser.parse_args()

    generateNetCDF(args.outfile, args.ndim, args.core_radius, args.gamma, args.x_center,
                   args.y_center, args.u_advection, args.v_advection, args.uz0, args.theoretical_model)
