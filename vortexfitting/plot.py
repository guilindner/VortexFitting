#!/usr/bin/env/ python
"""
Plotting routines and image generation
"""

import numpy as np
import matplotlib.pyplot as plt

from . import tools
from . import fitting


def plot_fields(vfield, detection_field):
    """
    Plot fields: display the (u,v,w) fields and the vorticity field.

    :param vfield: contains spatial mesh and velocity components
    :type vfield: class VelocityField()
    :param detection_field: detection field (vorticity ...)
    :type detection_field: 2D array of float
    :returns: popup
    :rtype: image
    """
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)  # , sharex='col', sharey='row')
    ax1.imshow(vfield.u_velocity_matrix, cmap='seismic', origin="lower")
    ax1.set_title('Velocity u (velocity_s)')

    ax2.imshow(vfield.v_velocity_matrix, cmap='seismic', origin="lower")
    ax2.set_title('Velocity v (velocity_n)')

    try:
        ax3.imshow(vfield.w_velocity_matrix, cmap='seismic', origin="lower")
        ax3.set_title('Velocity w (velocity_z)')
    except:
        print('No w velocity')
    ax4.set_title('Vorticity')
    ax4.imshow(detection_field, origin="lower", cmap='seismic')
    plt.tight_layout()

    plt.show()


def plot_detect(vortices_counterclockwise, vortices_clockwise, detection_field, *args):
    """
    Plot detect: display the location and rotation of the vortices

    :param vortices_counterclockwise: vortices spinning counterclockwise
    :type vortices_counterclockwise: list of vortices
    :param vortices_clockwise: vortices spinning clockwise
    :type vortices_clockwise: list of vortices
    :param detection_field: detection field (vorticity ...)
    :type detection_field: 2D array of float
    :returns: popup
    :rtype: image
    """
    plt.subplot()
    if (args[0] == True):
        detection_field = detection_field.T  # transpose the detection field
        plt.scatter(vortices_counterclockwise[0], vortices_counterclockwise[1], edgecolor='G', facecolor='G',
                    label='left')
        plt.scatter(vortices_clockwise[0], vortices_clockwise[1], edgecolor='Y', facecolor='Y', label='right')
    else:
        plt.scatter(vortices_counterclockwise[1], vortices_counterclockwise[0], edgecolor='G', facecolor='G',
                    label='left')
        plt.scatter(vortices_clockwise[1], vortices_clockwise[0], edgecolor='Y', facecolor='Y', label='right')

    plt.title('Detected possible vortices')
    # plt.contourf(field, cmap="Greys_r")

    plt.imshow(detection_field, origin='lower', cmap="Greys_r")
    plt.xlabel('x')
    plt.ylabel('y')
    # plt.legend()
    # plt.imshow(field, cmap="Greys_r",origin="lower")
    plt.tight_layout()

    plt.show()


def plot_quiver(x_index, y_index, u_data, v_data, detection_field):
    """
    Plot quiver: display a specific (x,y) location with vector fields.

    :param x_index: contains spatial mesh (x direction)
    :type x_index: array of float
    :param y_index: contains spatial mesh (y direction)
    :type y_index: array of float
    :param u_data: contains velocity data (u component)
    :type u_data: 2D array of float
    :param v_data: contains velocity data (v component)
    :type v_data: 2D array of float
    :param detection_field: detection field (vorticity ...)
    :type detection_field: 2D array of float
    :returns: popup
    :rtype: image
    """

    plt.figure()
    # plt.title('Velocity vectors centered at max swirling strength')
    plt.contourf(detection_field,
                 extent=[x_index[0][0], x_index[0][-1], y_index[0][0], y_index[-1][0]])
    s = 1  # sampling factor, can be modified
    plt.quiver(x_index[::s, ::s], y_index[::s, ::s], u_data[::s, ::s], v_data[::s, ::s])
    plt.show()


def plot_fit(x_index, y_index, u_data, v_data, u_model, v_model, 
             xc, yc, core_radius, gamma, u_advection, v_advection, correlation_value,
             cpt_vortex, subtract_advection_field, output_dir, time_step):
    """
    Plot fit

    :param x_index: contains spatial mesh (x direction)
    :type x_index: array of float
    :param y_index: contains spatial mesh (y direction)
    :type y_index: array of float
    :param u_data: contains velocity data (u component)
    :type u_data: 2D array of float
    :param v_data: contains velocity data (v component)
    :type v_data: 2D array of float
    :param u_model: contains velocity data (u component)
    :type u_model: 2D array of float
    :param v_model: contains velocity data (v component)
    :type v_model: 2D array of float
    :param xc: x coordinate of the vortex center
    :type xc: float
    :param yc: y coordinate of the vortex center
    :type yc: float
    :param core_radius: dimension of the vortex core radius
    :type core_radius: float
    :param gamma: circulation of the vortex
    :type gamma: float
    :param u_advection: contains velocity data (u component)
    :type u_advection: 2D array of float
    :param v_advection: contains velocity data (v component)
    :type v_advection: 2D array of float
    :param cpt_vortex: current n° of the vortex
    :type cpt_vortex: int
    :param subtract_advection_field: if True, the advection field (u_advection, v_advection) is subtracted
    :type subtract_advection_field:  bool
    :param output_dir: directory where the results are written
    :type output_dir: str
    :param correlation_value: correlation between the vortex and a Lamb-Oseen model
    :type correlation_value: float
    :param time_step: current time_step
    :type time_step: int
    :returns: image file
    :rtype: png
    """
    plt.figure()
    s = 1  # sampling factor
    if (x_index.size > 400):
        s = 2
    plt.quiver(x_index[::s, ::s], y_index[::s, ::s], u_data[::s, ::s], v_data[::s, ::s],
               color='r', label='data')
    plt.quiver(x_index[::s, ::s], y_index[::s, ::s], u_model[::s, ::s], v_model[::s, ::s],
               color='b', label='model', alpha=0.5)
    circle1 = plt.Circle((xc, yc), core_radius, color='k', alpha=0.05)
    plt.gca().add_artist(circle1)
    plt.gca().scatter([xc], [yc], marker='+', color='k', s=100)
    plt.legend()
    plt.grid()
    plt.gca().set_aspect('equal', adjustable='box')
    #    plt.axes().set_aspect('equal') #deprecated
    plt.xlabel('x')
    plt.ylabel('y')
    
    plt.title(r'r=%s $\Gamma$=%s u=%s v=%s C=%s' % (
    round(core_radius, 2), round(gamma, 2), round(u_advection, 2), round(v_advection, 2), round(correlation_value, 2)))
    if subtract_advection_field == False:
        plt.savefig(output_dir + '/vortex%i_%s.png' % (cpt_vortex, 'initial_vfield'), format='png')
    else:
        plt.savefig(output_dir + '/vortex%i_%s.png' % (cpt_vortex, 'advection_field_subtracted'), format='png')
    plt.close('all')


def plot_accepted(vfield, vortices_list, detection_field, output_dir, time_step):
    """
    Plot accepted: display the accepted vortices, with respect to the different criterion
    (correlation threshold, box size ...)

    :param vfield: contains spatial mesh and velocity components
    :type vfield: class VelocityField()
    :param vortices_list: contains all the detected vortices
    :type vortices_list: list
    :param detection_field: detection field (vorticity ...)
    :type detection_field: 2D array of float
    :param output_dir: directory where the results are written
    :type output_dir: str
    :param time_step: current time_step
    :type time_step: int
    :returns: popup
    :rtype: image
    """
    plt.figure(1)
    plt.contourf(vfield.x_coordinate_matrix, vfield.y_coordinate_matrix, detection_field, origin='lower', cmap="bone")
    plt.xlabel('x')
    plt.ylabel('y')
    dx = vfield.x_coordinate_step
    dy = vfield.y_coordinate_step
    plt.figure(2)
    plt.imshow(detection_field, origin='lower', cmap="bone")
    for i, line in enumerate(vortices_list):
        if vortices_list[i][1] > 0:
            orient = 'Y'
        else:
            orient = 'Y'
        plt.figure(1)
        circle1 = plt.Circle((line[2], line[3]), line[0],
                             edgecolor='yellow', facecolor='none', gid='vortex%i' % i)
        plt.gca().add_artist(circle1)
        plt.figure(2)
        circle1 = plt.Circle((line[2] / dx, line[3] / dy), line[0] / np.hypot(dx, dy),
                             edgecolor='yellow', facecolor='none', gid='vortex%i' % i)
        plt.gca().add_artist(circle1)

    ##Comparing data
    # fileIn = open('../data/dazin.dat', 'r')
    # for line in fileIn:
    # xComp = int(float(line.split()[1]))
    # yComp = int(float(line.split()[2]))
    # gammaComp = float(line.split()[3])
    # rComp = float(line.split()[4])
    # if gammaComp > 0:
    # orient = 'R'
    # else:
    # orient = 'R'
    # circle2=plt.Circle((xComp,yComp),rComp,edgecolor=orient,facecolor='none')
    # plt.gca().add_artist(circle2)

    # plt.legend()
    plt.figure(1)
    plt.tight_layout()
    plt.savefig(output_dir + '/accepted_{:01d}.svg'.format(time_step), format='svg')
    create_links(output_dir + '/accepted_{:01d}.svg'.format(time_step),vortices_list,output_dir,time_step)
    plt.figure(2)
    plt.savefig(output_dir + '/meshed_{:01d}.svg'.format(time_step),
                format='svg')  # use it to verify your coordinate system if needed !
    # plt.savefig(output_dir+'/tk_{:01d}.png'.format(time_step), format='png', transparent=True)
    
    # plt.show()


def plot_vortex(vfield, vortices_list, output_dir, time_step):
    """
    Plot vortex: plot a vortex and its corresponding vortex model

    :param vfield: contains spatial mesh and velocity components
    :type vfield: class VelocityField()
    :param vortices_list: contains all the detected vortices
    :type vortices_list: list
    :param output_dir: directory where the results are written
    :type output_dir: str
    :param time_step: current time_step
    :type time_step: int
    :returns: file
    :rtype: image
    """
    for cpt_vortex, line in enumerate(vortices_list):
        print('r: {:.3f}'.format(line[0]),
              'gamma: {:.2f}'.format(line[1]),
              'x: {:.2f}'.format(line[2]),
              'y: {:.2f}'.format(line[3]),
              'dist: {:.2f}'.format(line[6]),
              'corr: {:.2f}'.format(line[7]),
              'vt: {:.2f}'.format(line[8]))
        dx = vfield.x_coordinate_step
        dy = vfield.y_coordinate_step
        x_index, y_index, u_data, v_data = tools.window(vfield, round(line[2] / dx, 0), round(line[3] / dy, 0), line[6])
        u_model, v_model = fitting.velocity_model(line[0], line[1], line[2], line[3], line[4], line[5], x_index,
                                                  y_index)
        correlation_value = fitting.correlation_coef(u_data, v_data, u_model, v_model)
        plot_fit(x_index, y_index, u_data, v_data, u_model, v_model, line[2], line[3],
                 line[0], line[1], line[4], line[5], correlation_value, cpt_vortex, False, output_dir, time_step)
        correlation_value = fitting.correlation_coef(u_data - line[4], v_data - line[5], u_model - line[4], v_model - line[5])
        plot_fit(x_index, y_index, u_data - line[4], v_data - line[5], u_model - line[4], v_model - line[5], line[2],
                 line[3],
                 line[0], line[1], line[4], line[5], correlation_value, cpt_vortex, True, output_dir, time_step)


def create_links(path, vortices_list, output_dir, time_step):
    """
    create links: add some links bewteen the accepted.svg file and the detected vortices

    :param path: path of the accepted.svg file
    :type path: str
    :param vortices_list: contains all the detected vortices
    :type vortices_list: list
    :param output_dir: directory where the results are written
    :type output_dir: str
    :param time_step: current time_step
    :type time_step: int
    :returns: file
    :rtype: image
    """

    file_in = open(path, "r")
    file_out = open(output_dir + "/linked_{:01d}.svg".format(time_step), "w")
    i = 0
    vortex_found = False
    for line in file_in:
        if "</g>" in line:
            if vortex_found == True:
                file_out.write(line)
                file_out.write('   </a>\n')
                vortex_found = False
            else:
                file_out.write(line)
        elif "vortex" in line:
            file_out.write('   <a href="vortex%i_advection_field_subtracted.png">\n' % i)
            file_out.write(line)
            file_out.write('   <title>Vortex %i: r = %s gamma = %s</title>\n' % (
            i, round(vortices_list[i][3], 1), round(vortices_list[i][2], 1)))
            i = i + 1
            vortex_found = True
        else:
            file_out.write(line)
