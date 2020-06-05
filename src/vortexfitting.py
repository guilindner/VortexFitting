#!/usr/bin/env/ python3
"""vortex detection tool, by Guilherme Lindner, 2017-04\n
This program load NetCDF files from DNS simulations or PIV experiments, or Tecplot format.
It detects the vortices and apply a fitting to them.
"""
import sys
import argparse
import time
import numpy as np

from classes import VelocityField
import tools
import fitting
import plot
import schemes
import detection
import output

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Optional app description',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', dest='input_filename',
                        default='../data/test_dataHIT.nc',
                        help='input file', metavar='FILE')

    parser.add_argument('-o', '--output', dest='output_directory',
                        default='../results',
                        help='if you want to specify an output directory', metavar='DIRECTORY')

    parser.add_argument('-s', '--scheme', dest='scheme', type=int, default=22,
                        help='Scheme for differencing\n'
                             '2 = second order \n'
                             '22 = least-square filter (default)\n'
                             '4 = fourth order')

    parser.add_argument('-d', '--detect', dest='detection_method',
                        default='swirling',
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
                        default='/', 
                        help='Subtract a mean field from your instant field')

    parser.add_argument('-p', '--plot', dest='plot_method',
                        default='fit',
                        help='Plot on screen:\n'
                             'fit    = Detection and fitting, saves images (default)\n'
                             'detect = Possible vortices (no fitting)\n'
                             'fields = Velocity fields and vorticity\n')
                             
    parser.add_argument('-xy', '--xy', nargs=2, dest='xy_location', default=[0,0],
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
                        default='dns', help='Type of the file (default: dns)')

    parser.add_argument('-ct', '--corrthreshold', dest='correlation_threshold',
                        default=0.75, type=float, 
                        help='Correlation threshold (default: 0.75).' 
			     'if the vortex is too big, its better to decrease this value')

#    parser.add_argument('-v', '--verbose', dest='verbose',
#                        default=False, type=bool,
#                        help='Displays info or hides it. (default: True) ')

    args = parser.parse_args()

    start = time.time()

    #create file 'vortices.dat' to store output vortices data        
    output.create(args.output_directory) 

    #---- LOAD DATA ----#

    if (args.last < args.first):
        args.last= args.first

    for time_step in range(args.first,args.last+1,args.step):
        
        print('\nOpening file: ',args.input_filename.format(time_step), ', file type: ', args.file_type)
        if (args.mean_filename != '/'):
            print('Opening mean field: ', args.mean_filename)
        vfield = VelocityField(args.input_filename,time_step,args.mean_filename,args.file_type)

        #---- DIFFERENCE APPROXIMATION ----#
        lap = time.time()
        if args.scheme == 4:
            vfield.derivative = schemes.fourth_order_diff(vfield)
        elif args.scheme == 2:
            vfield.derivative = schemes.second_order_diff(vfield)
        elif args.scheme == 22:
            vfield.derivative = schemes.least_square_diff(vfield)
        else:
            print('No scheme', args.scheme, 'found. Exitting!')
            sys.exit()
        #print(round(time.time() - lap,3), 'seconds')
    
        #---- VORTICITY ----#
    
        vorticity = vfield.derivative['dvdx'] - vfield.derivative['dudy']
    
        #---- METHOD FOR DETECTION OF VORTICES ----#
        lap = time.time()
        if args.detection_method == 'Q':
            swirling = detection.q_criterion(vfield)
        elif args.detection_method == 'swirling':
            swirling = detection.calc_swirling(vfield)
        elif args.detection_method == 'delta':
            swirling = detection.delta_criterion(vfield)
        #print(round(time.time() - lap,3), 'seconds')
    
        if vfield.normalization_flag == True:
            swirling = tools.normalize(swirling,vfield.normalization_direction) #normalization
    
        #---- PEAK DETECTION ----#
        print('Threshold=',args.detection_threshold,', box size=',args.box_size)
    
        peaks = tools.find_peaks(swirling, args.detection_threshold, args.box_size)
    
        print('Vortices found: ',len(peaks[0]))
        #---- PEAKS DIRECTION OF ROTATION ----#
        vortices_counterclockwise, vortices_clockwise = tools.direction_rotation(vorticity,peaks)
    
        #---- MODEL FITTING ----#
        vortices = list()
        if (args.plot_method == 'fit' ) and (args.xy_location == [0, 0]):
            vortices = fitting.get_vortices(vfield,peaks,vorticity,args.rmax,args.correlation_threshold)
            print('---- Accepted vortices ----')
            print(len(vortices))
        else:
            print('No fitting')
    
        #---- PLOTTING OPTIONS ----#
        if args.xy_location != [0,0]:
            x_location = int(args.xy_location[0])
            y_location = int(args.xy_location[1])
            swirlingw = swirling[y_location-10:y_location+10,x_location-10:x_location+10]
            x_index, y_index, u_data, v_data = tools.window(vfield,x_location,y_location,10)
            plot.plot_quiver(x_index, y_index, u_data, v_data, swirlingw)
        if args.plot_method == 'detect':
            plot.plot_detect(vortices_counterclockwise,vortices_clockwise,swirling,args.flip_axis)
        if args.plot_method == 'fields':
            plot.plot_fields(vfield,vorticity)
        if args.plot_method == 'fit':
            plot.plot_accepted(vfield,vortices,swirling,args.output_directory,time_step)
            plot.plot_vortex(vfield,vortices,args.output_directory,time_step)
            output.write(vortices,args.output_directory,time_step)
