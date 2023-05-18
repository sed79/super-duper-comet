#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 15:58:16 2023

@author: sed79

Makes grid images of enhanced vs nonenhanced gas maps for each 
C2, NH2, and CN separately.

Enhanced maps are division by azimuthal median using https://www.psi.edu/research/cometimen on median coadded or cropped maps.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import simple_norm, ZScaleInterval, LinearStretch, ImageNormalize, ManualInterval, PercentileInterval
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve

def show_grid_CN():
    direc0 = '../data/MedianCoaddedMaps/CN_ALL/' #contains 400x400 maps of CN
    direc1 = '../data/EnhancedMaps/CN_ALL/' #only containing enhanced maps for 11-26, 12-05, 12-06, 12-21, 12-23, 12-29, 12-31, 03-16, 03-19
    for i, filename0 in enumerate(sorted(os.listdir(direc0))):
        for j, filename1 in enumerate(sorted(os.listdir(direc1))):
            if filename1[3:13] == filename0[3:13]:
                print(filename0[3:13])
                file0 = direc0+filename0
                file1 = direc1+filename1
                date = filename0[3:13]
                axindex = j % 3 
                print(axindex)
                if axindex == 0:
                    fig, ax = plt.subplots(3, 3, figsize=(10.5,9))
                #---unenhanced---
                hdulist = fits.open(file0,ignore_missing_end=True)
                image = hdulist[0].data
            
                norm = ImageNormalize(data=image, interval=ZScaleInterval(), stretch=LinearStretch())
                im = ax[axindex][0].pcolorfast(image, norm=norm, cmap='Greys_r')
                ax[axindex][0].set(xlim=(80,320),ylim=(80,320))
                ax[axindex][0].scatter(201, 201, marker='+', s=80, color='red')
                fig.colorbar(im, ax=ax[axindex][0], fraction=0.03, shrink=1, location='right', pad=0.02)
                
                # Species annotation
                props = dict(boxstyle='square', facecolor='k', alpha=0.7, pad=0.2, edgecolor='k')
                ax[axindex][0].text(95, 295, 'CN', color='white',fontsize=14, fontweight='bold',bbox=props)
                
                # Date annotation
                ax[axindex][0].text(200, 295, date, color='white',fontsize=14, fontweight='bold',bbox=props)
                
                # Type annotation
                ax[axindex][0].text(204, 270, 'Unenhanced', color='white',fontsize=12, fontweight='bold',bbox=props)
                
                # Get average angle values for the date
                table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
                names = table['new_arcfile_name']
                rows = table[names.str.contains(date)]
                velo = rows['velocityPA']
                sun = rows['sunTargetPA']
                avvelo = 360-np.mean(velo) #degrees clockwise from north
                avsun = 360-np.mean(sun) #degrees clockwise from north
                
                # antisunward arrow
                start = (90,90)
                l = 50
                dx = l * np.sin(avsun* np.pi / 180)
                dy = l * np.cos(avsun* np.pi / 180) 
                ax[axindex][0].arrow(start[0], start[1], dx, dy,fc='xkcd:fire engine red', ec='k',head_width=7,head_length=7,width=4, length_includes_head=True)
                # -S annotation
                txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
                props = dict(boxstyle='square', facecolor='black', alpha=0.7, pad=0.2, edgecolor='k')
                ax[axindex][0].text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=12, fontweight='bold', bbox=props) #annotation for arrow
                
                # antivelocity arrow
                start = (90,90)
                l = 50
                dx = l * np.sin(avvelo* np.pi / 180)
                dy = l * np.cos(avvelo* np.pi / 180) 
                ax[axindex][0].arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=7,head_length=7,width=4,length_includes_head=True)
                # -v annotation
                txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
                ax[axindex][0].text(txtpos[0], txtpos[1], '-V', color='white',fontsize=10, fontweight='bold', bbox = props) #annotation for arrow
            
                #Scale annotation
                dates = ['2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
                scale_date = [31.0997061, 32.9682025, 34.209806268718864, 34.328958743007014, 35.47341127722916, 35.532706858016816, 35.59324266213401, 35.57226744800711, 33.07054419420382, 32.84287518811179, 31.377380152136006, 30.253009538172453, 29.87299630886403, 27.75608220392316, 27.38499010989345] #from table of 10kkm apertures
                width = scale_date[dates.index(date)]
                centrepoint = 260
                ax[axindex][0].arrow(centrepoint-width/2, 95, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=4,length_includes_head=True)
                ax[axindex][0].text(220, 105, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
            
                #---enhanced---
                hdulist = fits.open(file1,ignore_missing_end=True)
                image = hdulist[0].data
                copyim = np.copy(image)
                for x in range(copyim.shape[0]):
                    for y in range(copyim.shape[1]):
                        if np.sqrt((x-200)**2+(y-200)**2) > 120:
                            copyim[x,y] = 'Nan'
                if date == '2019-12-21':
                    norm = simple_norm(copyim, 'linear', min_percent=0, max_percent=99.95)
                else:
                    norm = simple_norm(copyim, 'linear', min_percent=0, max_percent=99.95)
                im = ax[axindex][1].pcolorfast(copyim, norm=norm, cmap='Greys_r')
                ax[axindex][1].set(xlim=(80,320),ylim=(80,320))
                ax[axindex][1].scatter(201, 201, marker='+', s=80, color='red')
                fig.colorbar(im, ax=ax[axindex][1], fraction=0.03, shrink=1, location='right', pad=0.02)
                
                # Type annotation
                ax[axindex][1].text(100,100, 'Enhanced', color='white',fontsize=12, fontweight='bold',bbox=props)
            
                #---smoothed enhanced---
                numlevels = 10
                stddev = 5
                
                kernel = Gaussian2DKernel(x_stddev=stddev)
                convimage = convolve(image, kernel)
                copyconv = np.copy(convimage)
                for x in range(copyconv.shape[0]):
                    for y in range(copyconv.shape[1]):
                        if np.sqrt((x-200)**2+(y-200)**2) > 120:
                            copyconv[x,y] = 'Nan'
            
                norm = simple_norm(copyconv, 'linear', min_percent=0, max_percent=100)
                im = ax[axindex][2].pcolorfast(copyconv, norm=norm, cmap='Greys_r')
                cont = ax[axindex][2].contour(copyconv, levels=np.linspace(convimage.min(), convimage.max(), numlevels), colors='lime', origin='lower')
                ax[axindex][2].set(xlim=(80,320),ylim=(80,320))
                ax[axindex][2].scatter(201, 201, marker='+', s=80, color='red', zorder=5)
                fig.colorbar(im, ax=ax[axindex][2], fraction=0.03, shrink=1, location='right', pad=0.02)
                #ax[2].clabel(cont, cont.levels, inline=True, fmt=fmt, fontsize=6, colors='k')
                
                # Type annotation
                ax[axindex][2].text(100,100, 'Enhanced', color='white',fontsize=12, fontweight='bold',bbox=props)
                
                ax[axindex][0].axis('off')
                ax[axindex][1].axis('off')
                ax[axindex][2].axis('off')
                plt.tight_layout()
                
                if axindex == 2 or date == '2020-03-19':
                    print('saving')
                    plt.savefig('../output_images/Smooth_CN_grid_'+date+'.pdf')
    
    return

def show_grid_C2():
    direc0 = '../data/MedianCoaddedMaps/C2_ALL/'
    direc1 = '../data/EnhancedMaps/C2_ALL/'
    for i, filename0 in enumerate(sorted(os.listdir(direc0))):
        for j, filename1 in enumerate(sorted(os.listdir(direc1))):
            if filename1[3:13] == filename0[3:13]:
                print(filename0[3:13])
                file0 = direc0+filename0
                file1 = direc1+filename1
                date = filename0[3:13]
                axindex = j % 3 
                print(axindex)
                if axindex == 0:
                    fig, ax = plt.subplots(3, 3, figsize=(10.5,9))
                #---unenhanced---
                hdulist = fits.open(file0,ignore_missing_end=True)
                image = hdulist[0].data
            
                norm = ImageNormalize(data=image, interval=ZScaleInterval(), stretch=LinearStretch())
                im = ax[axindex][0].pcolorfast(image,norm=norm, cmap='Greys_r')
                ax[axindex][0].set(xlim=(80,320),ylim=(80,320))
                ax[axindex][0].scatter(201, 201, marker='+', s=80, color='red')
                fig.colorbar(im, ax=ax[axindex][0], fraction=0.03, shrink=1, location='right', pad=0.02)
                
                # Species annotation
                props = dict(boxstyle='square', facecolor='k', alpha=0.7, pad=0.2, edgecolor='k')
                ax[axindex][0].text(95, 295, 'C$_2$', color='white',fontsize=14, fontweight='bold',bbox=props)
                
                # Date annotation
                ax[axindex][0].text(200, 295, date, color='white',fontsize=14, fontweight='bold',bbox=props)
                
                # Type annotation
                ax[axindex][0].text(204, 270, 'Unenhanced', color='white',fontsize=12, fontweight='bold',bbox=props)
                
                # Get average angle values for the date
                table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
                names = table['new_arcfile_name']
                rows = table[names.str.contains(date)]
                velo = rows['velocityPA']
                sun = rows['sunTargetPA']
                avvelo = 360-np.mean(velo) #degrees clockwise from north
                avsun = 360-np.mean(sun) #degrees clockwise from north
                
                # antisunward arrow
                start = (90,90)
                l = 50
                dx = l * np.sin(avsun* np.pi / 180)
                dy = l * np.cos(avsun* np.pi / 180) 
                ax[axindex][0].arrow(start[0], start[1], dx, dy,fc='xkcd:fire engine red', ec='k',head_width=7,head_length=7,width=4, length_includes_head=True)
                # -S annotation
                txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
                props = dict(boxstyle='square', facecolor='black', alpha=0.7, pad=0.2, edgecolor='k')
                ax[axindex][0].text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=12, fontweight='bold', bbox=props) #annotation for arrow
                
                # antivelocity arrow
                start = (90,90)
                l = 50
                dx = l * np.sin(avvelo* np.pi / 180)
                dy = l * np.cos(avvelo* np.pi / 180) 
                ax[axindex][0].arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=7,head_length=7,width=4,length_includes_head=True)
                # -v annotation
                txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
                ax[axindex][0].text(txtpos[0], txtpos[1], '-V', color='white',fontsize=10, fontweight='bold', bbox = props) #annotation for arrow
            
                #Scale annotation
                dates = ['2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
                scale_date = [31.0997061, 32.9682025, 34.209806268718864, 34.328958743007014, 35.47341127722916, 35.532706858016816, 35.59324266213401, 35.57226744800711, 33.07054419420382, 32.84287518811179, 31.377380152136006, 30.253009538172453, 29.87299630886403, 27.75608220392316, 27.38499010989345] #from table of 10kkm apertures
                width = scale_date[dates.index(date)]
                centrepoint = 260
                ax[axindex][0].arrow(centrepoint-width/2, 95, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=4,length_includes_head=True)
                ax[axindex][0].text(220, 105, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
            
                #---enhanced---
                hdulist = fits.open(file1,ignore_missing_end=True)
                image = hdulist[0].data
                copyim = np.copy(image)
                for x in range(copyim.shape[0]):
                    for y in range(copyim.shape[1]):
                        if np.sqrt((x-200)**2+(y-200)**2) > 120:
                            copyim[x,y] = 'Nan'
                if date == '2019-12-21':
                    norm = simple_norm(copyim, 'linear', min_percent=0, max_percent=99.90)
                    #copyim[np.where(copyim > image.max()*0.9990)] = 'Nan' line didn't work
                else:
                    norm = simple_norm(copyim, 'linear', min_percent=0, max_percent=99.95)
                im = ax[axindex][1].pcolorfast(copyim, norm=norm, cmap='Greys_r')
                ax[axindex][1].set(xlim=(80,320),ylim=(80,320))
                ax[axindex][1].scatter(201, 201, marker='+', s=80, color='red')
                fig.colorbar(im, ax=ax[axindex][1], fraction=0.03, shrink=1, location='right', pad=0.02)
                
                # Type annotation
                ax[axindex][1].text(100,100, 'Enhanced', color='white',fontsize=12, fontweight='bold',bbox=props)
            
                #---smoothed enhanced---
                numlevels = 10
                stddev = 5
                
                kernel = Gaussian2DKernel(x_stddev=stddev)
                convimage = convolve(image, kernel)
                copyconv = np.copy(convimage)
                for x in range(copyconv.shape[0]):
                    for y in range(copyconv.shape[1]):
                        if np.sqrt((x-200)**2+(y-200)**2) > 120:
                            copyconv[x,y] = 'Nan'
                                 
                norm = simple_norm(copyconv, 'linear', min_percent=0, max_percent=100)
                im = ax[axindex][2].pcolorfast(copyconv, norm=norm, cmap='Greys_r')
                cont = ax[axindex][2].contour(copyconv, levels=np.linspace(convimage.min(), convimage.max(), numlevels), colors='lime', origin='lower')
                ax[axindex][2].set(xlim=(80,320),ylim=(80,320))
                ax[axindex][2].scatter(201, 201, marker='+', s=80, color='red', zorder=5)
                fig.colorbar(im, ax=ax[axindex][2], fraction=0.03, shrink=1, location='right', pad=0.02)
                #ax[2].clabel(cont, cont.levels, inline=True, fmt=fmt, fontsize=6, colors='k')
                
                # Type annotation
                ax[axindex][2].text(100,100, 'Enhanced', color='white',fontsize=12, fontweight='bold',bbox=props)
                
                ax[axindex][0].axis('off')
                ax[axindex][1].axis('off')
                ax[axindex][2].axis('off')
                plt.tight_layout()
                
                if axindex == 2 or date == '2020-03-19':
                    print('saving')
                    plt.savefig('../output_images/Smooth_C2_grid_'+date+'.pdf')
    
    return

def show_grid_NH2():
    direc0 = '../data/MedianCoaddedMaps/NH2_ALL/'
    direc1 = '..data/EnhancedMaps/NH2_ALL/'
    for i, filename0 in enumerate(sorted(os.listdir(direc0))):
        for j, filename1 in enumerate(sorted(os.listdir(direc1))):
            if filename1[4:14] == filename0[4:14]:
                print(filename0[4:14])
                file0 = direc0+filename0
                file1 = direc1+filename1
                date = filename0[4:14]
                axindex = j % 3 
                print(axindex)
                if axindex == 0:
                    fig, ax = plt.subplots(3, 3, figsize=(10.5,9))
                #---unenhanced---
                hdulist = fits.open(file0,ignore_missing_end=True)
                image = hdulist[0].data
            
                norm = ImageNormalize(data=image, interval=ZScaleInterval(), stretch=LinearStretch())
                im = ax[axindex][0].pcolorfast(image, norm=norm, cmap='Greys_r')
                ax[axindex][0].set(xlim=(80,320),ylim=(80,320))
                ax[axindex][0].scatter(201, 201, marker='+', s=80, color='red')
                fig.colorbar(im, ax=ax[axindex][0], fraction=0.03, shrink=1, location='right', pad=0.02)
                
                # Species annotation
                props = dict(boxstyle='square', facecolor='k', alpha=0.7, pad=0.2, edgecolor='k')
                ax[axindex][0].text(95, 295, 'NH$_2$', color='white',fontsize=14, fontweight='bold',bbox=props)
                
                # Date annotation
                ax[axindex][0].text(200, 295, date, color='white',fontsize=14, fontweight='bold',bbox=props)
                
                # Type annotation
                ax[axindex][0].text(204, 270, 'Unenhanced', color='white',fontsize=12, fontweight='bold',bbox=props)
                
                # Get average angle values for the date
                table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
                names = table['new_arcfile_name']
                rows = table[names.str.contains(date)]
                velo = rows['velocityPA']
                sun = rows['sunTargetPA']
                avvelo = 360-np.mean(velo) #degrees clockwise from north
                avsun = 360-np.mean(sun) #degrees clockwise from north
                
                # antisunward arrow
                start = (90,90)
                l = 50
                dx = l * np.sin(avsun* np.pi / 180)
                dy = l * np.cos(avsun* np.pi / 180) 
                ax[axindex][0].arrow(start[0], start[1], dx, dy,fc='xkcd:fire engine red', ec='k',head_width=7,head_length=7,width=4, length_includes_head=True)
                # -S annotation
                txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
                props = dict(boxstyle='square', facecolor='black', alpha=0.7, pad=0.2, edgecolor='k')
                ax[axindex][0].text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=12, fontweight='bold', bbox=props) #annotation for arrow
                
                # antivelocity arrow
                start = (90,90)
                l = 50
                dx = l * np.sin(avvelo* np.pi / 180)
                dy = l * np.cos(avvelo* np.pi / 180) 
                ax[axindex][0].arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=7,head_length=7,width=4,length_includes_head=True)
                # -v annotation
                txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
                ax[axindex][0].text(txtpos[0], txtpos[1], '-V', color='white',fontsize=10, fontweight='bold', bbox = props) #annotation for arrow
            
                #Scale annotation
                dates = ['2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
                scale_date = [31.0997061, 32.9682025, 34.209806268718864, 34.328958743007014, 35.47341127722916, 35.532706858016816, 35.59324266213401, 35.57226744800711, 33.07054419420382, 32.84287518811179, 31.377380152136006, 30.253009538172453, 29.87299630886403, 27.75608220392316, 27.38499010989345] #from table of 10kkm apertures
                width = scale_date[dates.index(date)]
                centrepoint = 260
                ax[axindex][0].arrow(centrepoint-width/2, 95, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=4,length_includes_head=True)
                ax[axindex][0].text(220, 105, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
            
                #---enhanced---
                hdulist = fits.open(file1,ignore_missing_end=True)
                image = hdulist[0].data
                copyim = np.copy(image)
                for x in range(copyim.shape[0]):
                    for y in range(copyim.shape[1]):
                        if np.sqrt((x-200)**2+(y-200)**2) > 120:
                            copyim[x,y] = 'Nan'
                if date == '2019-12-21':
                    norm = simple_norm(copyim, 'linear', min_percent=0, max_percent=99.90)
                    #copyim[np.where(copyim > image.max()*0.9990)] = 'Nan' line didn't work
                else:
                    norm = simple_norm(copyim, 'linear', min_percent=0, max_percent=99.95)
                im = ax[axindex][1].pcolorfast(copyim, norm=norm, cmap='Greys_r')
                ax[axindex][1].set(xlim=(80,320),ylim=(80,320))
                ax[axindex][1].scatter(201, 201, marker='+', s=80, color='red')
                fig.colorbar(im, ax=ax[axindex][1], fraction=0.03, shrink=1, location='right', pad=0.02)
                
                # Type annotation
                ax[axindex][1].text(100,100, 'Enhanced', color='white',fontsize=12, fontweight='bold',bbox=props)
            
                #---smoothed enhanced---
                numlevels = 10
                stddev = 5
                
                kernel = Gaussian2DKernel(x_stddev=stddev)
                convimage = convolve(image, kernel)
                copyconv = np.copy(convimage)
                for x in range(copyconv.shape[0]):
                    for y in range(copyconv.shape[1]):
                        if np.sqrt((x-200)**2+(y-200)**2) > 120:
                            copyconv[x,y] = 'Nan'
                                 
                norm = simple_norm(copyconv, 'linear', min_percent=0, max_percent=100)
                im = ax[axindex][2].pcolorfast(copyconv, norm=norm, cmap='Greys_r')
                cont = ax[axindex][2].contour(copyconv, levels=np.linspace(convimage.min(), convimage.max(), numlevels), colors='lime', origin='lower')
                ax[axindex][2].set(xlim=(80,320),ylim=(80,320))
                ax[axindex][2].scatter(201, 201, marker='+', s=80, color='red', zorder=5)
                fig.colorbar(im, ax=ax[axindex][2], fraction=0.03, shrink=1, location='right', pad=0.02)
                #ax[2].clabel(cont, cont.levels, inline=True, fmt=fmt, fontsize=6, colors='k')
                
                # Type annotation
                ax[axindex][2].text(100,100, 'Enhanced', color='white',fontsize=12, fontweight='bold',bbox=props)
                
                ax[axindex][0].axis('off')
                ax[axindex][1].axis('off')
                ax[axindex][2].axis('off')
                plt.tight_layout()
                
                if axindex == 2 or date == '2020-03-19':
                    print('saving')
                    plt.savefig('../output_images/Smooth_NH2_grid_'+date+'.pdf')
    
    return

def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"

show_grid_C2()

