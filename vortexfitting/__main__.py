#!/usr/bin/env/ python3

import sys
import os.path
import argparse

from vortexfitting import fitting
from vortexfitting import schemes
from vortexfitting import detection
from vortexfitting import output
from vortexfitting import classes


def main():
    parser = argparse.ArgumentParser(description='Optional app description',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', dest='input_filename',
                        default='data/example_dataHIT.nc', type=str,
                        help='Input file', metavar='FILE')

    parser.add_argument('-o', '--output', dest='output_directory',
                        default='results', type=str,
                        help='To specify an output directory', metavar='DIRECTORY')

    parser.add_argument('-s', '--scheme', dest='scheme',
                        default=22, type=int,
                        help='Scheme for differencing\n'
                             '2 = second order \n'
                             '22 = least-square filter (default)\n'
                             '4 = fourth order')

    parser.add_argument('-d', '--detect', dest='detection_method',
                        default='swirling', type=str,
                        help='Detection method:\n'
                             'Q = Q criterion\n'
                             'delta = delta criterion\n'
                             'swirling = 2D Swirling Strength (default)')

    parser.add_argument('-t', '--threshold', dest='detection_threshold',
                        default=0.0, type=float,
                        help='Threshold for detection (default=0.0)')

    parser.add_argument('-b', '--boxsize', dest='box_size',
                        default=6, type=int,
                        help='Box size for the detection (default=6)')

    parser.add_argument('-f', '--flip', dest='flip_axis',
                        default=False, type=bool,
                        help='Flip X and Y axis for plotting:\n'
                             '0 = False (default)\n'
                             '1 = True')

    parser.add_argument('-mf', '--meanfilename', dest='mean_filename',
                        default='/', type=str,
                        help='Subtract a mean field from your instant field')

    parser.add_argument('-p', '--plot', dest='plot_method',
                        default='fit', type=str,
                        help='Plot on screen:\n'
                             'fit    = Detection and fitting, saves images (default)\n'
                             'detect = Possible vortices (no fitting)\n'
                             'fields = Velocity fields and vorticity\n')

    parser.add_argument('-xy', '--xy', nargs=2, dest='xy_location',
                        default=[0, 0],
                        help='specify a location to see the data. ex: -xy 80 60')

    parser.add_argument('-first', '--first', dest='first',
                        default=0, type=int,
                        help='Index of first field (default: 0)')

    parser.add_argument('-last', '--last', dest='last',
                        default=0, type=int,
                        help='Index of last field (default: 0)')

    parser.add_argument('-step', '--step', dest='step',
                        default=1, type=int,
                        help='Step between 2 fields (default: 1)')

    parser.add_argument('-rmax', '--rmax', dest='rmax',
                        default=10, type=float,
                        help='Initial guess on the vortex radius (default: 10)')

    parser.add_argument('-ft', '--filetype', dest='file_type',
                        default='dns', type=str,
                        help='Type of the file (default: dns)\n'
                             'piv_netcdf: read a netCDF format\n'
                             'dns: read a netCDF format\n'
                             'piv_tecplot: read a tecplot format\n'
                             'openfoam: read an openfoam format\n')

    parser.add_argument('-ct', '--corrthreshold', dest='correlation_threshold',
                        default=0.75, type=float,
                        help='Correlation threshold (default: 0.75).\n'
                             'If the vortex is too big, its better to decrease this value')

    parser.add_argument('-of', '--outputformat', dest='output_format',
                        default='png', type=str,
                        help='Format of the output file (default: png).\n'
                             'Can be png, pdf, jpg ...')
    #    parser.add_argument('-v', '--verbose', dest='verbose',
    #                        default=False, type=bool,
    #                        help='Displays info or hides it. (default: True) ')

    args = parser.parse_args()
    # start = time.time()

    # ---- LOAD DATA ----#

    if args.last < args.first:
        args.last = args.first  # last should be at least equal to first

    for time_step in range(args.first, args.last + 1, args.step):

        if not os.path.exists(args.input_filename.format(time_step)):
            print('The input file does not exist. Exiting.')
            sys.exit()

        print('\nOpening file: ', args.input_filename.format(time_step), ', file type: ', args.file_type)
        if args.mean_filename != '/':
            print('Opening mean field: ', args.mean_filename)

        vfield = classes.VelocityField(args.input_filename, time_step, args.mean_filename, args.file_type)

        # ---- DIFFERENCE APPROXIMATION ----#
        # lap = time.time()
        if args.scheme == 4:
            vfield.derivative = schemes.fourth_order_diff(vfield)
        elif args.scheme == 2:
            vfield.derivative = schemes.second_order_diff(vfield)
        elif args.scheme == 22:
            vfield.derivative = schemes.least_square_diff(vfield)
        else:
            print('No scheme', args.scheme, 'found. Exiting!')
            sys.exit()
        # print(round(time.time() - lap,3), 'seconds')

        # ---- VORTICITY ----#

        vorticity = vfield.derivative['dvdx'] - vfield.derivative['dudy']

        # ---- METHOD FOR DETECTION OF VORTICES ----#
        # lap = time.time()

        detection_field = []
        if args.detection_method == 'Q':
            detection_field = detection.calc_q_criterion(vfield)
        elif args.detection_method == 'swirling':
            detection_field = detection.calc_swirling(vfield)
        elif args.detection_method == 'delta':
            detection_field = detection.calc_delta_criterion(vfield)

        # print(round(time.time() - lap,3), 'seconds')

        if vfield.normalization_flag:
            print('Normalization for ', vfield.normalization_direction, ' direction')
            detection_field = fitting.normalize(detection_field, vfield.normalization_direction)  # normalization

        # ---- PEAK DETECTION ----#
        print('Threshold=', args.detection_threshold, ', box size=', args.box_size)

        peaks = fitting.find_peaks(detection_field, args.detection_threshold, args.box_size)

        print('Vortices found: ', len(peaks[0]))
        # ---- PEAKS DIRECTION OF ROTATION ----#
        vortices_counterclockwise, vortices_clockwise = fitting.direction_rotation(vorticity, peaks)

        # ---- MODEL FITTING ----#
        vortices = list()
        if (args.plot_method == 'fit') and (args.xy_location == [0, 0]):
            vortices = fitting.get_vortices(vfield, peaks, vorticity, args.rmax, args.correlation_threshold)
            print('---- Accepted vortices ----')
            print(len(vortices))
        else:
            print('No fitting')

        # ---- PLOTTING OPTIONS ----#
        if args.xy_location != [0, 0]:
            x_location = int(args.xy_location[0])
            y_location = int(args.xy_location[1])
            detection_field_window = detection_field[y_location - 10:y_location + 10, x_location - 10:x_location + 10]
            x_index, y_index, u_data, v_data = fitting.window(vfield, x_location, y_location, 10)
            fitting.plot_quiver(x_index, y_index, u_data, v_data, detection_field_window)
        if args.plot_method == 'detect':
            fitting.plot_detect(vortices_counterclockwise, vortices_clockwise, detection_field, args.flip_axis)
        if args.plot_method == 'fields':
            fitting.plot_fields(vfield, vorticity)
        if args.plot_method == 'fit':
            # create file 'vortices.dat' to store output vortices data        
            output.create(args.output_directory, args)
            fitting.plot_accepted(vfield, vortices, detection_field, args.output_directory, time_step,
                                  args.output_format)
            fitting.plot_vortex(vfield, vortices, args.output_directory, time_step, args.output_format)
            output.write(vortices, args.output_directory, time_step)
            # output.write_field(args.output_directory+'/detection_field.dat', args.detection_method,
            #                    vfield,detection_field)
