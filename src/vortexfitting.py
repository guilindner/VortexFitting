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
                        default='../data/test_data.nc',
                        help='input NetCDF file', metavar='FILE')
                        
    parser.add_argument('-o', '--output', dest='outfilename',
                        help='output NetCDF file', metavar='FILE')
    
    parser.add_argument('-s', '--scheme', dest='scheme', type=int, default=2,
                        help='Scheme for differencing\n'
                             '2 = second order\n'
                             '4 = fourth order')
    
    parser.add_argument('-T', '--time', dest='timestep', type=int,
                        default=0,
                        help='Timestep/Sample desired')
                        
    parser.add_argument('-d', '--detect', dest='detect',
                        default='swirling',
                        help='Detection method:\n'
                             'Q = Q criterion\n'
                             'swirling = 2D Swirling Strength')
    
    parser.add_argument('-t', '--threshold', dest='threshold',
                        default=1.5, type=float,
                        help='Threshold for detection, integer')

    parser.add_argument('-b', '--boxsize', dest='boxsize',
                        default=6, type=int,
                        help='Box size for the detection')
    
    parser.add_argument('-f', '--flip', dest='flip',
                        default=False, type=bool,
                        help='Flip X and Y axis for plotting, 0 = False, 1 = True')
    
    parser.add_argument('-p', '--plot', dest='plot_x',
                        default='',
                        help='Plot on screen:\n'
                             'detect = Vortices position\n'
                             'fields = Velocity fields\n'
                             'quiver = Vector on specific position')
    
    args = parser.parse_args()
    
    start = time.time()
    #---- LOAD DATA ----#
    print("Opening file:",args.infilename)

    #print("Sample target: (todo)", args.timestep)
    
    a = VelocityField(args.infilename,args.timestep)
    print("Samples:", a.samples)

    #---- DIFFERENCE APPROXIMATION ----# 
    lap = time.time()
    if args.scheme == 4:
        a.derivative = schemes.fourth_order_diff(a)
    elif args.scheme == 2:
        a.derivative = schemes.second_order_diff(a)
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
    #print(round(time.time() - lap,3), 'seconds')

    if a.norm == True:
        swirling = tools.normalize(swirling,a.normdir) #normalization

    #---- PEAK DETECTION ----#
    print("threshold=",args.threshold,"box size=",args.boxsize)

    peaks = tools.find_peaks(swirling, args.threshold, args.boxsize)

    print("Vortices found:",len(peaks[0]))

    #---- PEAKS DIRECTION OF ROTATION ----#
    dirL, dirR = tools.direction_rotation(vorticity,peaks)

    #---- MODEL FITTING ----# SEE IN PLOT
    vortices = list()
    
    for i in range(len(peaks[0])):
            xCenter = peaks[0][i]
            yCenter = peaks[1][i]
            if (244 > xCenter > 10) and (244 > yCenter > 10):
                gamma = vorticity[xCenter,yCenter]
                coreR, gamma, corr, dist, fxCenter, fyCenter, u_conv, v_conv = fitting.full_fit(a, xCenter, yCenter, gamma)
                #print(a.dx[xCenter],fxCenter,'|',a.dy[yCenter],fyCenter)
                if (corr > 0.75):
                    #2vortices.append([xCenter,yCenter, gamma, coreR,corr,dist])
                    vortices.append([xCenter,yCenter, gamma, coreR,corr,dist,fxCenter,fyCenter,u_conv,v_conv]) #not fitted to plot the center!  
    print('---- Accepted vortices ----')
    print(len(vortices))
    #print('xCenter, yCenter, gamma, core Radius, correlation, mesh distance')
    #for vortex in vortices:
        #print(vortex)

    #---- SAVING OUTPUT FILE ----#
    if args.outfilename == None:
        pass
    else:
        print("saving file",args.outfilename)
    
  
    #---- PLOTTING OPTIONS ----#
    if args.plot_x == 'detect':
        plot.plot_detection(dirL,dirR,swirling,args.flip)
    elif args.plot_x == 'fields':
        plot.plot_fields(a,vorticity)
    elif args.plot_x == 'quiverRuim':
        dist = 10
        for i in range(len(peaks[0])):
            xCenter = peaks[0][i]
            yCenter = peaks[1][i]
            X, Y, Uw, Vw = tools.window(a,xCenter,yCenter,dist)
            swirlingw = swirling[xCenter-dist:xCenter+dist,yCenter-dist:yCenter+dist] #reuse window function?
            if (xCenter > dist) and (yCenter > dist):
                print('x1:',xCenter,'x2:',yCenter, 'swirl:',peaks[2][i])
                plot.plot_quiver(X, Y, Uw, Vw, swirlingw)
    elif args.plot_x == 'quiver':
        for i in range(len(vortices)):
            xCenter = vortices[i][0]
            yCenter = vortices[i][1]
            gamma = vortices[i][2]
            coreR = vortices[i][3]
            corr = vortices[i][4]
            dist = vortices[i][5]
            u_conv = vortices[i][6]
            v_conv = vortices[i][7]
            swirlingw = swirling[xCenter-dist:xCenter+dist,yCenter-dist:yCenter+dist]
            X, Y, Uw, Vw = tools.window(a,xCenter,yCenter,dist)
            uMod, vMod = fitting.velocity_model(u_conv, v_conv, X, Y,fxCenter,fyCenter, gamma, coreR)
            plot.plot_quiver(X, Y, Uw, Vw, swirlingw)
                
    elif args.plot_x == 'fit':
        for i in range(len(vortices)):
            xCenter = vortices[i][0]
            yCenter = vortices[i][1]
            gamma = vortices[i][2]
            coreR = vortices[i][3]
            corr = vortices[i][4]
            dist = vortices[i][5]
            fxCenter = vortices[i][6]
            fyCenter = vortices[i][7]
            u_conv = vortices[i][8]
            v_conv = vortices[i][9]
            dx = a.dx[xCenter+1]-a.dx[xCenter]
            dy = a.dy[yCenter+1]-a.dx[yCenter]
            fshiftxCenter = (a.dx[xCenter] -fxCenter)/dx
            fshiftyCenter = (a.dy[yCenter] -fyCenter)/dy
            shiftxCenter = int(fshiftxCenter)
            shiftyCenter = int(fshiftyCenter)
            print(xCenter, shiftxCenter)
            print(a.dx[xCenter], fxCenter)
            xCenter = xCenter - shiftxCenter
            yCenter = yCenter - shiftyCenter

            print('xC:',xCenter,'yC:',yCenter, 'vort:',gamma, 'mesh',dist, 'corr',corr, 'coreR',coreR)
            X, Y, Uw, Vw = tools.window(a,xCenter,yCenter,dist)
            uMod, vMod = fitting.velocity_model(u_conv, v_conv, X, Y, fxCenter, fyCenter, gamma, coreR)
            corr = fitting.correlation_coef(Uw,Vw,uMod,vMod)
            plot.plot_corr(X, Y, Uw, Vw, uMod, vMod, coreR, corr)
    
    else:
        print('no plot')
