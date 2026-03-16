#!/usr/bin/env/ python3
"""
Create an output file for the detected vortices, with tecplot format
"""

import os
import numpy as np
from classes import VelocityField


def create(output_directory: str, args: list) -> None:
    """
    Create an output file

    :param output_directory: directory hosting the file vortices.dat
    :type output_directory: str
    :param args: args
    :type args: class parser
    :returns: file with time, radius, gamma, xcenter, ycenter, u_advection, v_advection, correlation, vtheta
    :rtype: file
    """
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    if not os.path.exists(output_directory + '/vortices.dat'):
        outfile = open(output_directory + '/vortices.dat', 'w')
        outfile.write("TITLE=\"Vortex characteristics evolution\"\n")
        outfile.write("Variables=\"time\",\"radius\",\"gamma\",\"xcenter\",\"ycenter\","
                      "\"u_advection\",\"v_advection\",\"correlation\",\"vtheta\",\"uz0\"\n")
        outfile.write("DATASETAUXDATA Detection_method=\"{}\"\n".format(args.detection_method))
        outfile.write("DATASETAUXDATA Scheme=\"{}\"\n".format(args.scheme))
        outfile.write("DATASETAUXDATA Box_size=\"{}\"\n".format(args.box_size))
        outfile.write("DATASETAUXDATA Detection_threshold=\"{}\"\n".format(args.detection_threshold))
        outfile.write("DATASETAUXDATA Rmax=\"{}\"\n".format(args.rmax))
        outfile.write("DATASETAUXDATA Correlation_threshold=\"{}\"\n".format(args.correlation_threshold))
        outfile.write("DATASETAUXDATA Vortex_Model=\"{}\"\n".format(args.theoretical_model))
        outfile.write("DATASETAUXDATA Mean_file=\"{}\"\n".format(args.mean_filename))
        outfile.write("DATASETAUXDATA File_type=\"{}\"\n".format(args.file_type))
        outfile.write("ZONE T=\"0\", SOLUTIONTIME=0\n")
        outfile.close()


def write(vortices: list, output_directory: str, time_step: int) -> None:
    """
    Update an output file

    :param vortices: list of the detected vortices
    :param output_directory: directory hosting the file vortices.dat
    :param time_step: time of the current velocity field
    :type vortices: list 
    :type output_directory: str
    :type time_step: int
    :returns: empty
    :rtype: empty
    """

    outfile = open(output_directory + '/vortices.dat', 'a')

    # outfile.write("time r0 gamma xc yc ua va corr utheta uz0\n")
    for i, line in enumerate(vortices):
        outfile.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n".format(time_step, line[0], line[1], line[2], line[3],
                                                                         line[4], line[5], line[7], line[8], line[9]))
    outfile.close()


def write_field(output_file: str, detection_method: str, vfield: VelocityField, detection_field: np.ndarray) -> None:
    """
    Write a detection field file

    :param output_file: directory hosting the file vortices.dat
    :param detection_method: writes the selected detection method (Q, Delta, swirling strength)
    :param vfield: full size velocity field
    :param detection_field: full size detection field
    :type output_file: str
    :type detection_method: str
    :type vfield: ndarray
    :type detection_field: ndarray
    :returns: file
    :rtype: file
    """
    if not os.path.exists(output_file):
        os.makedirs(output_file)
    outfile = open(output_file, 'w')
    outfile.write("TITLE=\"Detection field\"\n")
    outfile.write("Variables=\"X\",\"Y\",\"{}\"\n".format(detection_method))
    outfile.write(
        "ZONE T=\"0\", I={:d}, J={:d}, SOLUTIONTIME=0\n".format(vfield.x_coordinate_size, vfield.y_coordinate_size))
    for j in np.arange(0, vfield.y_coordinate_size, 1):
        for i in np.arange(0, vfield.x_coordinate_size, 1):
            outfile.write("{0} {1} {2}\n".format(str(vfield.x_coordinate_matrix[j]), str(vfield.y_coordinate_matrix[i]),
                                                 detection_field[i, j]))
    outfile.write("{0} {1} {2}\n".format(0, 0, 0))
    outfile.close()
