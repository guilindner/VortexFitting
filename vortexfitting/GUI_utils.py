#!/usr/bin/env/ python3
"""
Some functions for the Graphical User Interface
"""

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import fitting
from classes import VelocityField


def gui_plot_fit(
        ax: matplotlib.axes.Axes,
        canvas: matplotlib.backends.backend_agg.FigureCanvasAgg,
        x_index: np.ndarray,
        y_index: np.ndarray,
        u_data: np.ndarray,
        v_data: np.ndarray,
        u_model: np.ndarray,
        v_model: np.ndarray,
        xc: float,
        yc: float,
        core_radius: float,
        gamma: float,
        u_advection: np.ndarray,
        v_advection: np.ndarray,
        correlation_value: float,
        flip: bool
) -> None:
    """
    Quiver plot of data velocity and model velocity in GUI

    :param ax: Matplotlib AxesSubplot for the current plot
    :param canvas: FigureCanvas associated to the axis
    :param x_index: spatial grid (x direction)
    :param y_index: spatial grid (y direction)
    :param u_data: data velocity component (u)
    :param v_data: data velocity component (v)
    :param u_model: model velocity component (u)
    :param v_model: model velocity component (v)
    :param xc: x  coordinate of the vortex center
    :param yc: y coordinate of the vortex center
    :param core_radius: vortex radius
    :param gamma: vortex circulation
    :param u_advection: advection velocity component (u)
    :param v_advection: advection velocity component (v)
    :param correlation_value: corrélation entre le vortex et un modèle Lamb-Oseen.
    :param flip: invert x and y axes if True.  
    :type ax: figure
    :type canvas: canvas
    :type x_index: float matrix
    :type y_index: float matrix
    :type u_data: float matrix
    :type v_data: float matrix
    :type u_model: float matrix
    :type v_model: float matrix
    :type xc: float
    :type yc: float
    :type core_radius: float
    :type gamma: float
    :type u_advection: float matrix
    :type v_advection: float matrix
    :type correlation_value: float
    :type flip: boolean

    :returns: empty
    :rtype: empty           
    """
    ax.clear()

    s = 1  # Facteur d'échantillonnage
    if x_index.size > 400:
        s = 2

    if flip:
        print("flip value", flip)
        ax.quiver(y_index[::s, ::s], x_index[::s, ::s], np.transpose(u_data[::s, ::s]),
                  np.transpose(v_data[::s, ::s]), color='r', label='data')
        ax.quiver(y_index[::s, ::s], x_index[::s, ::s], np.transpose(u_model[::s, ::s]),
                  np.transpose(v_model[::s, ::s]), color='b', label='model', alpha=0.5)
        circle = plt.Circle((yc, xc), core_radius, color='k', alpha=0.05)
        ax.add_artist(circle)
        ax.scatter([yc], [xc], marker='+', color='k', s=100)
        ax.set_xlabel('y')
        ax.set_ylabel('x')
    else:
        ax.quiver(x_index[::s, ::s], y_index[::s, ::s], u_data[::s, ::s], v_data[::s, ::s],
                  color='r', label='data')
        ax.quiver(x_index[::s, ::s], y_index[::s, ::s], u_model[::s, ::s], v_model[::s, ::s],
                  color='b', label='model', alpha=0.5)
        circle = plt.Circle((xc, yc), core_radius, color='k', alpha=0.05)
        ax.add_artist(circle)
        ax.scatter([xc], [yc], marker='+', color='k', s=100)
        ax.set_xlabel('x')
        ax.set_ylabel('y')

    ax.legend()
    ax.grid()
    ax.set_aspect('equal', adjustable='box')

    ax.set_title(r'r=%s $\Gamma$=%s u=%s v=%s C=%s' % (
        round(core_radius, 2), round(gamma, 2), np.round(u_advection, 2), np.round(v_advection, 2),
        round(correlation_value, 2)))

    canvas.draw()


def gui_plot_accepted(
        ax: matplotlib.axes.Axes,
        canvas: matplotlib.backends.backend_agg.FigureCanvasAgg,
        vfield: VelocityField,
        vortices_list: list,
        detection_field: np.ndarray,
        flip: bool) -> None:
    """
    Plot accepted vortices in GUI: display the accepted vortices with respect to the different criteria.

    :param ax: Matplotlib AxesSubplot for the current plot
    :type ax: figure
    :param canvas: FigureCanvas associated to the axis
    :type canvas: canva
    :param vfield: spatial grid and velocity components
    :param vortices_list: all the detected vortices
    :param detection_field: detection field (vorticity ...)
    :param flip: invert x and y axes if True.
    :type vfield: class vfield
    :type vortices_list: list
    :type detection_field: np.ndarray
    :type flip: boolean

    :returns: empty
    :rtype: empty
    """
    ax.clear()

    if flip:
        ax.contourf(vfield.y_coordinate_matrix, vfield.x_coordinate_matrix, np.transpose(detection_field),
                    origin='lower', cmap="bone")
        ax.set_ylabel('x')
        ax.set_xlabel('y')

        # dy = vfield.x_coordinate_step
        # dx = vfield.y_coordinate_step
        for i, line in enumerate(vortices_list):
            if vortices_list[i][1] > 0:  # Gamma > 0: counterclockwise
                circle = plt.Circle((line[3], line[2]), line[0], edgecolor='green', facecolor='none')
            else:  # Gamma < 0: clockwise
                circle = plt.Circle((line[3], line[2]), line[0], edgecolor='yellow', facecolor='none')
            ax.add_artist(circle)
        
    else:
        ax.contourf(vfield.x_coordinate_matrix, vfield.y_coordinate_matrix, detection_field, origin='lower',
                    cmap="bone")
        ax.set_xlabel('x')
        ax.set_ylabel('y')

        # dx = vfield.x_coordinate_step
        # dy = vfield.y_coordinate_step

        for i, line in enumerate(vortices_list):
            if vortices_list[i][1] > 0:  # Gamma > 0: counterclockwise
                circle = plt.Circle((line[2], line[3]), line[0], edgecolor='green', facecolor='none')
            else:  # Gamma < 0: clockwise
                circle = plt.Circle((line[2], line[3]), line[0], edgecolor='yellow', facecolor='none')
            ax.add_artist(circle)

    canvas.draw()


def gui_plot_vortex(
        ax: matplotlib.axes.Axes,
        canvas: matplotlib.backends.backend_agg.FigureCanvasAgg,
        vfield: VelocityField,
        vortices: list,
        theoretical_model: str,
        flip: bool
) -> None:
    """
    Plot vortex: plot a vortex and its corresponding vortex model

    :param ax: Matplotlib AxesSubplot for the current plot
    :type ax: figure
    :param canvas: FigureCanvas associated to the axis
    :type canvas: canva
    :param vfield: contains spatial mesh and velocity components
    :type vfield: class VelocityField
    :param vortices: contains all the detected vortices
    :type vortices: list
    :param theoretical_model: chosen theoretical model
    :type theoretical_model: str
    :param flip: invert x and y axes if True.
    :type flip: boolean

    :returns: file
    :rtype: image
    """

    dx = vfield.x_coordinate_step
    dy = vfield.y_coordinate_step

    if theoretical_model == 'batchelor':
        x_index, y_index, u_data, v_data, w_data = fitting.window(vfield, round(vortices[2] / dx, 0),
                                                                  round(vortices[3] / dy, 0), vortices[6],
                                                                  theoretical_model)
        u_model, v_model, w_model = fitting.velocity_model(vortices[0], vortices[1], vortices[2], vortices[3],
                                                           vortices[4], vortices[5], x_index,
                                                           y_index, np.mean(w_data), theoretical_model)
    else:
        x_index, y_index, u_data, v_data = fitting.window(vfield, round(vortices[2] / dx, 0),
                                                          round(vortices[3] / dy, 0), vortices[6], theoretical_model)
        u_model, v_model = fitting.velocity_model(vortices[0], vortices[1], vortices[2], vortices[3], vortices[4],
                                                  vortices[5], x_index,
                                                  y_index, None, theoretical_model)
    correlation_value = fitting.correlation_coef(u_data, v_data, None, u_model, v_model, None)

    gui_plot_fit(ax, canvas, x_index, y_index, u_data, v_data, u_model, v_model,
                 vortices[2], vortices[3], vortices[0], vortices[1], vortices[4], vortices[5], correlation_value, flip)

    print('r: {:.3f}'.format(vortices[0]),
          'gamma: {:.2f}'.format(vortices[1]),
          'xc: {:.2f}'.format(vortices[2]),
          'yc: {:.2f}'.format(vortices[3]),
          'correlation: {:.2f}'.format(correlation_value),
          'utheta: {:.2f}'.format(vortices[8]),
          'uz0: {:.2f}'.format(vortices[9]))
