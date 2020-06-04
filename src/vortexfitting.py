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

    parser.add_argument('-i', '--input', dest='infilename',
                        default='../data/test_dataHIT.nc',
                        help='input file', metavar='FILE')

    parser.add_argument('-o', '--output', dest='output_dir',
                        default='../results',
                        help='if you want to specify an output directory', metavar='DIRECTORY')

    parser.add_argument('-s', '--scheme', dest='scheme', type=int, default=22,
                        help='Scheme for differencing\n'
                             '2 = second order \n'
                             '22 = least-square filter (default)\n'
                             '4 = fourth order')

#    parser.add_argument('-T', '--time', dest='timestep', type=int,
#                        default=0,
#                        help='Timestep/Sample/Z position desired (default=0)')

    parser.add_argument('-d', '--detect', dest='detect',
                        default='swirling',
                        help='Detection method:\n'
                             'Q = Q criterion\n'
                             'delta = delta criterion\n'
                             'swirling = 2D Swirling Strength (default)')

    parser.add_argument('-t', '--threshold', dest='threshold',
                        default=0.0, type=float,
                        help='Threshold for detection (default=0.0)')

    parser.add_argument('-b', '--boxsize', dest='boxsize',
                        default=6, type=int,
                        help='Box size for the detection (default=6)')

    parser.add_argument('-f', '--flip', dest='flip',
                        default=False, type=bool,
                        help='Flip X and Y axis for plotting:\n'
                              '0 = False (default)\n'
                              '1 = True')

    parser.add_argument('-mf', '--meanfilename', dest='meanfilename',
                        default='/', 
                        help='Subtract a mean field from your instant field')

    parser.add_argument('-p', '--plot', dest='plot_x',
                        default='fit',
                        help='Plot on screen:\n'
                             'fit    = Detection and fitting, saves images (default)\n'
                             'detect = Possible vortices (no fitting)\n'
                             'fields = Velocity fields and vorticity\n')
                             
    parser.add_argument('-xy', '--xy', nargs=2, dest='xy', default=[0,0],
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

    parser.add_argument('-ft', '--filetype', dest='filetype',
                        default='dns', help='Type of the file (default: dns)')

    parser.add_argument('-ct', '--corrthreshold', dest='corr_threshold',
                        default=0.75, help='Correlation threshold (default: 0.75).' 
					     'if the vortex is too big, its better to decrease this value')

#    parser.add_argument('-v', '--verbose', dest='verbose',
#                        default=False, type=bool,
#                        help='Displays info or hides it. (default: True) ')

    args = parser.parse_args()

    start = time.time()

    #create file 'vortices.dat' to store output vortices data        
    output.create(args.output_dir) 

    #---- LOAD DATA ----#

    if (args.last < args.first):
        args.last= args.first

    for time_step in range(args.first,args.last+1,args.step):
        
        print('\nOpening file: ',args.infilename.format(time_step), ', file type: ', args.filetype)
        if (args.meanfilename != '/'):
            print('Opening mean field: ', args.meanfilename)
        vfield = VelocityField(args.infilename,time_step,args.meanfilename,args.filetype)
        
        # print("Samples:", vfield.samples)
    
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
        if args.detect == 'Q':
            swirling = detection.q_criterion(vfield)
        elif args.detect == 'swirling':
            swirling = detection.calc_swirling(vfield)
        elif args.detect == 'delta':
            swirling = detection.delta_criterion(vfield)
        #print(round(time.time() - lap,3), 'seconds')
    
        if vfield.norm == True:
            swirling = tools.normalize(swirling,vfield.normdir) #normalization
    
        #---- PEAK DETECTION ----#
        print('Threshold=',args.threshold,', box size=',args.boxsize)
    
        peaks = tools.find_peaks(swirling, args.threshold, args.boxsize)
    
        print('Vortices found: ',len(peaks[0]))
        #---- PEAKS DIRECTION OF ROTATION ----#
        dirL, dirR = tools.direction_rotation(vorticity,peaks)
    
        #---- MODEL FITTING ----#
        vortices = list()
        if (args.plot_x == 'fit' ) and (args.xy == [0, 0]):
            vortices = fitting.get_vortices(vfield,peaks,vorticity,args.rmax,args.corr_threshold)
            print('---- Accepted vortices ----')
            print(len(vortices))
        else:
            print('No fitting')
    
        #---- PLOTTING OPTIONS ----#
        if args.xy != [0,0]:
            x = int(args.xy[0])
            y = int(args.xy[1])
            swirlingw = swirling[y-10:y+10,x-10:x+10]
            x_index, y_index, u_data, v_data = tools.window(vfield,x,y,10)
            plot.plot_quiver(x_index, y_index, u_data, v_data, swirlingw)
        if args.plot_x == 'detect':
            plot.plot_detect(dirL,dirR,swirling,args.flip)
        if args.plot_x == 'fields':
            plot.plot_fields(vfield,vorticity)
        if args.plot_x == 'fit':
            plot.plot_accepted(vfield,vortices,swirling,args.output_dir,time_step)
            plot.plot_vortex(vfield,vortices,args.output_dir,time_step)
            output.write(vortices,args.output_dir,time_step)
