#!/usr/bin/env/ python3
"""vortex detection tool, by Guilherme Lindner, 2017-04\n
This program load NetCDF files from DNS simulations  or PIV experiments
and detect the vortices and apply a fitting to them.
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Optional app description',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input', dest='infilename',
                        default='../data/test_dataHIT.nc',
                        help='input file', metavar='FILE')

    parser.add_argument('-ft', '--filetype', dest='filetype',
                        default='dns',
                        help='file type: piv, dns (default), tecplot', metavar='FILETYPE')

    parser.add_argument('-o', '--output', dest='output_dir',
                        default='../results',
                        help='if you want to specify an output directory', metavar='DIRECTORY')

    parser.add_argument('-s', '--scheme', dest='scheme', type=int, default=22,
                        help='Scheme for differencing\n'
                             '2 = second order (default)\n'
                             '22 = least-square filter\n'
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
                        help='Average file tu subtract')

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
                        help='First image to read')

    parser.add_argument('-last', '--last', dest='last',
                        default=0, type=int,
                        help='Last image to read')

    parser.add_argument('-step', '--step', dest='step',
                        default=1, type=int,
                        help='Step between 2 images')

    parser.add_argument('-rmax', '--rmax', dest='rmax',
                        default=0.0, type=float,
                        help='guess on the starting vortex radius')

    args = parser.parse_args()

    start = time.time()
    #---- LOAD DATA ----#

	#init output vortices file, with tecplot format
    #Change someday (function of format...)
    outfile = open(args.output_dir+'/vortices.dat','w')
    outfile.write("TITLE=\"Vortex characteristics evolution\"\n")
    outfile.write("Variables=\"time\",\"radius\",\"gamma\",\"xindex\",\"yindex\",\"uc\",\"vc\",\"dist\",\"corr\",\"vtheta\"\n")
    outfile.write("ZONE T=\"0\", SOLUTIONTIME=0\n")
    outfile.close()

    for time_step in range(args.first,args.last+1,args.step):
        print
        a = VelocityField(args.infilename,time_step,args.meanfilename,args.filetype)
        
        print("Samples:", a.samples)
    
        #---- DIFFERENCE APPROXIMATION ----#
        lap = time.time()
        if args.scheme == 4:
            a.derivative = schemes.fourth_order_diff(a)
        elif args.scheme == 2:
            a.derivative = schemes.second_order_diff(a)
        elif args.scheme == 22:
            a.derivative = schemes.least_square_diff(a)
        else:
            print('No scheme', args.scheme, 'found. Exitting!')
            sys.exit()
        #print(round(time.time() - lap,3), 'seconds')
    
        #---- VORTICITY ----#
    
        vorticity = a.derivative['dvdx'] - a.derivative['dudy']
    
        #---- METHOD FOR DETECTION OF VORTICES ----#
        lap = time.time()
        if args.detect == 'Q':
            swirling = detection.q_criterion(a)
        elif args.detect == 'swirling':
            swirling = detection.calc_swirling(a)
        elif args.detect == 'delta':
            swirling = detection.delta_criterion(a)
        #print(round(time.time() - lap,3), 'seconds')
    
        if a.norm == True:
            swirling = tools.normalize(swirling,a.normdir) #normalization
    
        #---- PEAK DETECTION ----#
        print("threshold=",args.threshold,"box size=",args.boxsize)
    
        peaks = tools.find_peaks(swirling, args.threshold, args.boxsize)
    
        print("Vortices found:",len(peaks[0]))
        #---- PEAKS DIRECTION OF ROTATION ----#
        dirL, dirR = tools.direction_rotation(vorticity,peaks)
    
        #---- MODEL FITTING ----#
        vortices = list()
        if (args.plot_x == 'fit' ) and (args.xy == [0, 0]):
            vortices = fitting.get_vortices(a,peaks,vorticity,args.rmax)
            print('---- Accepted vortices ----')
            print(len(vortices))
        else:
            print("No fitting")
    
        #---- PLOTTING OPTIONS ----#
        if args.xy != [0,0]:
            x = int(args.xy[0])
            y = int(args.xy[1])
            swirlingw = swirling[y-10:y+10,x-10:x+10]
            x_index, y_index, u_data, v_data = tools.window(a,x,y,10)
            plot.plot_quiver(x_index, y_index, u_data, v_data, swirlingw)
        if args.plot_x == 'detect':
            plot.plot_detect(dirL,dirR,swirling,args.flip)
        if args.plot_x == 'fields':
            plot.plot_fields(a,vorticity)
        if args.plot_x == 'fit':
            plot.plot_accepted(a,vortices,swirling,args.output_dir,time_step)
            plot.plot_vortex(a,vortices,args.output_dir,time_step)
