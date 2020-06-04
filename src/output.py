#!/usr/bin/env/ python3
"""
Create an output file for the detected vortices, with tecplot format
"""
import os


def create(output_directory):
    """
    Create an output file

    :param output_directory: directory hosting the file vortices.dat
    :type output_directory: str
    """
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    outfile = open(output_directory+'/vortices.dat','w')
    outfile.write("TITLE=\"Vortex characteristics evolution\"\n")
    outfile.write("Variables=\"time\",\"radius\",\"gamma\",\"xindex\",\"yindex\","
                  "\"uc\",\"vc\",\"dist\",\"corr\",\"vtheta\"\n")
    outfile.write("ZONE T=\"0\", SOLUTIONTIME=0\n")
    outfile.close()


def write(vortices, output_directory, time_step):
    """
    Update an output file

    :param vortices: list of the detected vortices
    :param output_directory: directory hosting the file vortices.dat
    :param time_step: time of the current velocity field
    :type vortices: list 
    :type output_directory: str
    :type time_step: int
    """

    outfile = open(output_directory+'/vortices.dat', 'a')
    for i, line in enumerate(vortices):
        outfile.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n".format(time_step, line[0], line[1], line[2], line[3],
                                                                         line[4], line[5], line[6], line[7], line[8]))



