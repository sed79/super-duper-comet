#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 14:50:39 2023

@author: sed79

Plot uncombined datacubes summed over the spectral axis next to
CN emission maps for comparison of where stars are relative to 2I.

Separate plot for each month of observations. Only uses certain exposures per night,
so only these can be in directory for grid positioning of plots to work.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.text as mtext
from astropy.io import fits
from astropy.visualization import simple_norm, ZScaleInterval, LinearStretch, ImageNormalize, ManualInterval, PercentileInterval, SqrtStretch, AsinhStretch, SinhStretch
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.visualization import make_lupton_rgb
from mpdaf.obj import Cube

            
            
def grid_allflux_march():
    directory =  '../data/LambdaSum/' #folder of lambda sum cubes (cubes summed over spectral axis)
    fig, ax = plt.subplots(2, 5, figsize=(18.75,7))
    j=-1
    for i, filename in enumerate(sorted(os.listdir(directory))):
        if filename[18:26] == '2020-03-':
            j+=1
            
            if j <= 3:
                a=0
                b=j
            elif 4<=j<=6:
                a=1
                b=j-2

            date = filename[18:28]
            time = filename[29:37]

            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            imageuncropped = hdulist[0].data
            imageonecrop = imageuncropped[100:340, 90:330]
            image = imageonecrop[55:165,60:180]
            
            # display map on
            norm = ImageNormalize(data=image, interval=PercentileInterval(99.5), stretch=AsinhStretch()) #ZScaleInterval
            im = ax[a,b].pcolorfast(image, norm=norm, cmap='Greys_r')
            
            # Species annotation
            props = dict(boxstyle='square', facecolor='k', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(8, 100, 'All Flux', color='white',fontsize=12, fontweight='bold',bbox=props)
            
            # Date annotation
            ax[a,b].text(80, 100, date, color='white',fontsize=12, fontweight='bold',bbox=props)
            ax[a,b].text(90, 92, time, color='white',fontsize=10, fontweight='bold',bbox=props)
            
            # Get average angle values for the date
            table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
            names = table['new_arcfile_name']
            rows = table[names.str.contains(date)]
            velo = rows['velocityPA']
            sun = rows['sunTargetPA']
            avvelo = 360-np.mean(velo) #degrees clockwise from north
            avsun = 360-np.mean(sun) #degrees clockwise from north
            
            # antisunward arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avsun* np.pi / 180)
            dy = l * np.cos(avsun* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='dimgrey', ec='k',head_width=7,head_length=7,width=2, length_includes_head=True)
            # -S annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            props = dict(boxstyle='square', facecolor='black', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=12, fontweight='bold', bbox=props) #annotation for arrow
            
            # antivelocity arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avvelo* np.pi / 180)
            dy = l * np.cos(avvelo* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=7,head_length=7,width=2,length_includes_head=True)
            # -v annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            ax[a,b].text(txtpos[0], txtpos[1], '-V', color='white',fontsize=10, fontweight='bold', bbox = props) #annotation for arrow
        
    
            #Scale annotation
            dates = ['2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
            scale_date = [31.0997061, 32.9682025, 34.209806268718864, 34.328958743007014, 35.47341127722916, 35.532706858016816, 35.59324266213401, 35.57226744800711, 33.07054419420382, 32.84287518811179, 31.377380152136006, 30.253009538172453, 29.87299630886403, 27.75608220392316, 27.38499010989345] #from table of 10kkm apertures
            width = scale_date[dates.index(date)]
            centrepoint = 95
            ax[a,b].arrow(centrepoint-width/2, 10, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=2,length_includes_head=True)
            ax[a,b].text(81, 14, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
            
            ax[a,b].axis('off')
            ax[a,b].set_aspect('equal')
    directory =  '../CroppedMaps/CN_ALL/' #CN emission maps
    j=-1
    for i, filename in enumerate(sorted(os.listdir(directory))):
        if filename[11:19] == '2020-03-':
            j+=1
            a = j
            b = 4

            date = filename[11:21]
            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            imageuncropped = hdulist[0].data
            image = imageuncropped[55:165,60:180]
            
            # display map on zscale
            norm = ImageNormalize(data=image, interval=ZScaleInterval(), stretch=LinearStretch())
            im = ax[a,b].pcolorfast(image, norm=norm, cmap='viridis')
            kernel = Gaussian2DKernel(x_stddev=9)
            convimage = convolve(image, kernel)
            ax[a,b].contour(convimage, levels=np.linspace(convimage.min(), convimage.max(), 12), colors='whitesmoke', alpha=0.7, origin='lower', linewidths=1, zorder=1)
            
            # Species annotation
            props = dict(boxstyle='square', facecolor='k', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(8, 100, 'CN', color='white',fontsize=12, fontweight='bold',bbox=props)
            
            # Date annotation
            ax[a,b].text(80, 100, date, color='white',fontsize=12, fontweight='bold',bbox=props)
            
            # Get average angle values for the date
            table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
            names = table['new_arcfile_name']
            rows = table[names.str.contains(date)]
            velo = rows['velocityPA']
            sun = rows['sunTargetPA']
            avvelo = 360-np.mean(velo) #degrees clockwise from north
            avsun = 360-np.mean(sun) #degrees clockwise from north
            
            # antisunward arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avsun* np.pi / 180)
            dy = l * np.cos(avsun* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='dimgrey', ec='k',head_width=7,head_length=7,width=2, length_includes_head=True)
            # -S annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            props = dict(boxstyle='square', facecolor='black', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=12, fontweight='bold', bbox=props) #annotation for arrow
            
            # antivelocity arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avvelo* np.pi / 180)
            dy = l * np.cos(avvelo* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=7,head_length=7,width=2,length_includes_head=True)
            # -v annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            ax[a,b].text(txtpos[0], txtpos[1], '-V', color='white',fontsize=10, fontweight='bold', bbox = props) #annotation for arrow
        
    
            #Scale annotation
            dates = ['2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
            scale_date = [31.0997061, 32.9682025, 34.209806268718864, 34.328958743007014, 35.47341127722916, 35.532706858016816, 35.59324266213401, 35.57226744800711, 33.07054419420382, 32.84287518811179, 31.377380152136006, 30.253009538172453, 29.87299630886403, 27.75608220392316, 27.38499010989345] #from table of 10kkm apertures
            width = scale_date[dates.index(date)]
            centrepoint = 95
            ax[a,b].arrow(centrepoint-width/2, 10, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=2,length_includes_head=True)
            ax[a,b].text(81, 14, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
            
            ax[a,b].axis('off')
            ax[a,b].set_aspect('equal')
            
    ax[1, 0].axis('off')
    ax[1, 1].axis('off')
    plt.tight_layout(pad=0, w_pad=0, h_pad=0)
    ax[1, 0].remove()
    ax[1, 1].remove()

    bbox = fig.bbox_inches.from_bounds(-0.2, 0, 19,7)
    plt.savefig('../output_images/Compare_allflux_CN_elongation_march.pdf', bbox_inches=bbox)
    
    return
    

def grid_allflux_feb():
    directory =  '../data/LambdaSum/'
    fig, ax = plt.subplots(5, 5, figsize=(18.75,17.5))
    j=-1
    for i, filename in enumerate(sorted(os.listdir(directory))):
        if filename[18:26] == '2020-02-':
            j+=1
            
            if j <= 3:
                a=0
                b=j
            elif j == 4:
                a=1
                b=j-1
            elif 5<=j<=8:
                a=2
                b=j-5
            elif 9<=j<=12:
                a=3
                b=j-9
            elif 13<=j<=15:
                a=4
                b=j-12

            date = filename[18:28]
            time = filename[29:37]

            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            imageuncropped = hdulist[0].data
            imageonecrop = imageuncropped[100:340, 90:330]
            image = imageonecrop[55:165,60:180]
            
            # display map
            norm = ImageNormalize(data=image, interval=PercentileInterval(99.5), stretch=AsinhStretch()) #ZScaleInterval
            im = ax[a,b].pcolorfast(image, norm=norm, cmap='Greys_r')
            
            # Species annotation
            props = dict(boxstyle='square', facecolor='k', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(8, 100, 'All Flux', color='white',fontsize=12, fontweight='bold',bbox=props)
            
            # Date annotation
            ax[a,b].text(80, 100, date, color='white',fontsize=12, fontweight='bold',bbox=props)
            ax[a,b].text(90, 92, time, color='white',fontsize=10, fontweight='bold',bbox=props)
            
            # Get average angle values for the date
            table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
            names = table['new_arcfile_name']
            rows = table[names.str.contains(date)]
            velo = rows['velocityPA']
            sun = rows['sunTargetPA']
            avvelo = 360-np.mean(velo) #degrees clockwise from north
            avsun = 360-np.mean(sun) #degrees clockwise from north
            
            # antisunward arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avsun* np.pi / 180)
            dy = l * np.cos(avsun* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='dimgrey', ec='k',head_width=7,head_length=7,width=2, length_includes_head=True)
            # -S annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            props = dict(boxstyle='square', facecolor='black', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=14, fontweight='bold', bbox=props) #annotation for arrow
            
            # antivelocity arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avvelo* np.pi / 180)
            dy = l * np.cos(avvelo* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=7,head_length=7,width=2,length_includes_head=True)
            # -v annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            ax[a,b].text(txtpos[0], txtpos[1], '-V', color='white',fontsize=12, fontweight='bold', bbox = props) #annotation for arrow
        
    
            #Scale annotation
            dates = ['2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
            scale_date = [31.0997061, 32.9682025, 34.209806268718864, 34.328958743007014, 35.47341127722916, 35.532706858016816, 35.59324266213401, 35.57226744800711, 33.07054419420382, 32.84287518811179, 31.377380152136006, 30.253009538172453, 29.87299630886403, 27.75608220392316, 27.38499010989345] #from table of 10kkm apertures
            width = scale_date[dates.index(date)]
            centrepoint = 95
            ax[a,b].arrow(centrepoint-width/2, 10, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=2,length_includes_head=True)
            ax[a,b].text(81, 14, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
            
            ax[a,b].axis('off')
            ax[a,b].set_aspect('equal')
       
    directory =  '../CroppedMaps/CN_ALL/' #CN emission maps
    j=-1
    for i, filename in enumerate(sorted(os.listdir(directory))):
        if filename[11:19] == '2020-02-':
            j+=1
            a = j
            b = 4

            date = filename[11:21]
            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            imageuncropped = hdulist[0].data
            image = imageuncropped[55:165,60:180]
            
            # display map on zscale
            norm = ImageNormalize(data=image, interval=ZScaleInterval(), stretch=LinearStretch())
            im = ax[a,b].pcolorfast(image, norm=norm, cmap='viridis')
            kernel = Gaussian2DKernel(x_stddev=9)
            convimage = convolve(image, kernel)
            ax[a,b].contour(convimage, levels=np.linspace(convimage.min(), convimage.max(), 12), colors='whitesmoke',  alpha=0.7, origin='lower', linewidths=1, zorder=1)
            
            # Species annotation
            props = dict(boxstyle='square', facecolor='k', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(8, 100, 'CN', color='white',fontsize=12, fontweight='bold',bbox=props)
            
            # Date annotation
            ax[a,b].text(80, 100, date, color='white',fontsize=12, fontweight='bold',bbox=props)
            
            # Get average angle values for the date
            table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
            names = table['new_arcfile_name']
            rows = table[names.str.contains(date)]
            velo = rows['velocityPA']
            sun = rows['sunTargetPA']
            avvelo = 360-np.mean(velo) #degrees clockwise from north
            avsun = 360-np.mean(sun) #degrees clockwise from north
            
            # antisunward arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avsun* np.pi / 180)
            dy = l * np.cos(avsun* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='dimgrey', ec='k',head_width=7,head_length=7,width=2, length_includes_head=True)
            # -S annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            props = dict(boxstyle='square', facecolor='black', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=14, fontweight='bold', bbox=props) #annotation for arrow
            
            # antivelocity arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avvelo* np.pi / 180)
            dy = l * np.cos(avvelo* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=7,head_length=7,width=2,length_includes_head=True)
            # -v annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            ax[a,b].text(txtpos[0], txtpos[1], '-V', color='white',fontsize=12, fontweight='bold', bbox = props) #annotation for arrow
        
    
            #Scale annotation
            dates = ['2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
            scale_date = [31.0997061, 32.9682025, 34.209806268718864, 34.328958743007014, 35.47341127722916, 35.532706858016816, 35.59324266213401, 35.57226744800711, 33.07054419420382, 32.84287518811179, 31.377380152136006, 30.253009538172453, 29.87299630886403, 27.75608220392316, 27.38499010989345] #from table of 10kkm apertures
            width = scale_date[dates.index(date)]
            centrepoint = 95
            ax[a,b].arrow(centrepoint-width/2, 10, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=2,length_includes_head=True)
            ax[a,b].text(81, 14, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
            
            ax[a,b].axis('off')
            ax[a,b].set_aspect('equal')
            
    ax[1, 0].axis('off')
    ax[1, 1].axis('off')
    ax[1, 2].axis('off')
    ax[4, 0].axis('off')
    plt.tight_layout(pad=0, w_pad=0, h_pad=0)
    ax[1, 0].remove()
    ax[1, 1].remove()
    ax[1, 2].remove()
    ax[4, 0].remove()
    bbox = fig.bbox_inches.from_bounds(-0.2, 0, 19,17.5)
    plt.savefig('../output_images/Compare_allflux_CN_elongation_feb.pdf', bbox_inches=bbox)
    
    return
            
            
def grid_allflux_dec():
    directory =  '../data/LambdaSum/'
    fig, ax = plt.subplots(6, 5, figsize=(18.75,21))
    j=-1
    for i, filename in enumerate(sorted(os.listdir(directory))):
        if filename[18:26] == '2019-12-':
            j+=1
            
            if j <= 3:
                a=0
                b=j
            elif 4<=j<=6:
                a=1
                b=j-3
            elif 7<=j<=8:
                a=2
                b=j-5
            elif 9<=j<=12:
                a=3
                b=j-9
            elif 13<=j<=16:
                a=4
                b=j-13
            elif 17<=j<=20:
                a=5
                b=j-17

            date = filename[18:28]
            time = filename[29:37]

            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            imageuncropped = hdulist[0].data
            imageonecrop = imageuncropped[100:340, 90:330]
            image = imageonecrop[55:165,60:180]
            
            # display map
            norm = ImageNormalize(data=image, interval=PercentileInterval(99.5), stretch=AsinhStretch()) #ZScaleInterval
            im = ax[a,b].pcolorfast(image, norm=norm, cmap='Greys_r')
            
            # Species annotation
            props = dict(boxstyle='square', facecolor='k', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(8, 100, 'All Flux', color='white',fontsize=12, fontweight='bold',bbox=props)
            
            # Date annotation
            ax[a,b].text(80, 100, date, color='white',fontsize=12, fontweight='bold',bbox=props)
            ax[a,b].text(90, 92, time, color='white',fontsize=10, fontweight='bold',bbox=props)
            
            # Get average angle values for the date
            table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
            names = table['new_arcfile_name']
            rows = table[names.str.contains(date)]
            velo = rows['velocityPA']
            sun = rows['sunTargetPA']
            avvelo = 360-np.mean(velo) #degrees clockwise from north
            avsun = 360-np.mean(sun) #degrees clockwise from north
            
            # antisunward arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avsun* np.pi / 180)
            dy = l * np.cos(avsun* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='dimgrey', ec='k',head_width=7,head_length=7,width=2, length_includes_head=True)
            # -S annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            props = dict(boxstyle='square', facecolor='black', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=14, fontweight='bold', bbox=props) #annotation for arrow
            
            # antivelocity arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avvelo* np.pi / 180)
            dy = l * np.cos(avvelo* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=7,head_length=7,width=2,length_includes_head=True)
            # -v annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            ax[a,b].text(txtpos[0], txtpos[1], '-V', color='white',fontsize=12, fontweight='bold', bbox = props) #annotation for arrow
        
    
            #Scale annotation
            dates = ['2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
            scale_date = [31.0997061, 32.9682025, 34.209806268718864, 34.328958743007014, 35.47341127722916, 35.532706858016816, 35.59324266213401, 35.57226744800711, 33.07054419420382, 32.84287518811179, 31.377380152136006, 30.253009538172453, 29.87299630886403, 27.75608220392316, 27.38499010989345] #from table of 10kkm apertures
            width = scale_date[dates.index(date)]
            centrepoint = 95
            ax[a,b].arrow(centrepoint-width/2, 10, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=2,length_includes_head=True)
            ax[a,b].text(81, 14, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
            
            ax[a,b].axis('off')
            ax[a,b].set_aspect('equal')

    directory =  '../CroppedMaps/CN_ALL/' #CN emission maps
    j=-1
    for i, filename in enumerate(sorted(os.listdir(directory))):
        if filename[11:19] == '2019-12-':
            j+=1
            a = j
            b = 4

            date = filename[11:21]
            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            imageuncropped = hdulist[0].data
            image = imageuncropped[55:165,60:180]
            
            # display map on zscale
            norm = ImageNormalize(data=image, interval=ZScaleInterval(), stretch=LinearStretch())
            im = ax[a,b].pcolorfast(image, norm=norm, cmap='viridis')
            kernel = Gaussian2DKernel(x_stddev=9)
            convimage = convolve(image, kernel)
            ax[a,b].contour(convimage, levels=np.linspace(convimage.min(), convimage.max(), 12), colors='whitesmoke', alpha=0.7, origin='lower', linewidths=1, zorder=1)
            
            # Species annotation
            props = dict(boxstyle='square', facecolor='k', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(8, 100, 'CN', color='white',fontsize=12, fontweight='bold',bbox=props)
            
            # Date annotation
            ax[a,b].text(80, 100, date, color='white',fontsize=12, fontweight='bold',bbox=props)
            
            # Get average angle values for the date
            table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
            names = table['new_arcfile_name']
            rows = table[names.str.contains(date)]
            velo = rows['velocityPA']
            sun = rows['sunTargetPA']
            avvelo = 360-np.mean(velo) #degrees clockwise from north
            avsun = 360-np.mean(sun) #degrees clockwise from north
            
            # antisunward arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avsun* np.pi / 180)
            dy = l * np.cos(avsun* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='dimgrey', ec='k',head_width=7,head_length=7,width=2, length_includes_head=True)
            # -S annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            props = dict(boxstyle='square', facecolor='black', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=14, fontweight='bold', bbox=props) #annotation for arrow
            
            # antivelocity arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avvelo* np.pi / 180)
            dy = l * np.cos(avvelo* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=7,head_length=7,width=2,length_includes_head=True)
            # -v annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            ax[a,b].text(txtpos[0], txtpos[1], '-V', color='white',fontsize=12, fontweight='bold', bbox = props) #annotation for arrow
        
    
            #Scale annotation
            dates = ['2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
            scale_date = [31.0997061, 32.9682025, 34.209806268718864, 34.328958743007014, 35.47341127722916, 35.532706858016816, 35.59324266213401, 35.57226744800711, 33.07054419420382, 32.84287518811179, 31.377380152136006, 30.253009538172453, 29.87299630886403, 27.75608220392316, 27.38499010989345] #from table of 10kkm apertures
            width = scale_date[dates.index(date)]
            centrepoint = 95
            ax[a,b].arrow(centrepoint-width/2, 10, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=2,length_includes_head=True)
            ax[a,b].text(81, 14, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
            
            ax[a,b].axis('off')
            ax[a,b].set_aspect('equal')
    
    ax[1, 0].axis('off')
    ax[2, 1].axis('off')
    ax[2, 0].axis('off')
    plt.tight_layout(pad=0, w_pad=0, h_pad=0)    
    ax[1, 0].remove()
    ax[2, 1].remove()
    ax[2, 0].remove()

    bbox = fig.bbox_inches.from_bounds(-0.2, 0, 19,21)
    plt.savefig('../output_images/Compare_allflux_CN_elongation_dec.pdf', bbox_inches=bbox)
    
    return
          
            
            
def grid_allflux_nov():
    directory =  '../data/LambdaSum/'
    fig, ax = plt.subplots(3, 5, figsize=(18.75,10.5))
    j=-1
    for i, filename in enumerate(sorted(os.listdir(directory))):
        if filename[18:26] == '2019-11-':
            j+=1
            
            if j <= 1:
                a=0
                b=j+2
            elif 2<=j<=3:
                a=1
                b=j
            elif 4<=j:
                a=2
                b=j-4

            date = filename[18:28]
            time = filename[29:37]

            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            imageuncropped = hdulist[0].data
            imageonecrop = imageuncropped[100:340, 90:330]
            image = imageonecrop[55:165,60:180]
            
            # display map on zscale
            
            norm = ImageNormalize(data=image, interval=PercentileInterval(99.5), stretch=AsinhStretch()) #ZScaleInterval
            im = ax[a,b].pcolorfast(image, norm=norm, cmap='Greys_r')
             
            # Centre point
            #plt.scatter(50.5, 45.5, s=100, c='red', marker='+')
            
            # Species annotation
            props = dict(boxstyle='square', facecolor='k', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(8, 100, 'All Flux', color='white',fontsize=12, fontweight='bold',bbox=props)
            
            # Date annotation
            ax[a,b].text(80, 100, date, color='white',fontsize=12, fontweight='bold',bbox=props)
            ax[a,b].text(90, 92, time, color='white',fontsize=10, fontweight='bold',bbox=props)
            
            # Get average angle values for the date
            table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
            names = table['new_arcfile_name']
            rows = table[names.str.contains(date)]
            velo = rows['velocityPA']
            sun = rows['sunTargetPA']
            avvelo = 360-np.mean(velo) #degrees clockwise from north
            avsun = 360-np.mean(sun) #degrees clockwise from north
            
            # antisunward arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avsun* np.pi / 180)
            dy = l * np.cos(avsun* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='dimgrey', ec='k',head_width=7,head_length=7,width=2, length_includes_head=True)
            # -S annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            props = dict(boxstyle='square', facecolor='black', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=12, fontweight='bold', bbox=props) #annotation for arrow
            
            # antivelocity arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avvelo* np.pi / 180)
            dy = l * np.cos(avvelo* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=7,head_length=7,width=2,length_includes_head=True)
            # -v annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            ax[a,b].text(txtpos[0], txtpos[1], '-V', color='white',fontsize=10, fontweight='bold', bbox = props) #annotation for arrow
        
    
            #Scale annotation
            dates = ['2019-11-14', '2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
            scale_date = [31.0997061, 31.0997061, 32.9682025, 34.209806268718864, 34.328958743007014, 35.47341127722916, 35.532706858016816, 35.59324266213401, 35.57226744800711, 33.07054419420382, 32.84287518811179, 31.377380152136006, 30.253009538172453, 29.87299630886403, 27.75608220392316, 27.38499010989345] #from table of 10kkm apertures
            width = scale_date[dates.index(date)]
            centrepoint = 95
            ax[a,b].arrow(centrepoint-width/2, 10, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=2,length_includes_head=True)
            ax[a,b].text(81, 14, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
            
            ax[a,b].axis('off')
            ax[a,b].set_aspect('equal')
            
      
       
    directory =  '../CroppedMaps/CN_ALL/' #CN emission maps
    j=-1
    for i, filename in enumerate(sorted(os.listdir(directory))):
        if filename[11:19] == '2019-11-':
            j+=1
            a = j
            b = 4

            date = filename[11:21]
            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            if date == '2019-11-14':
                imageuncropped = hdulist[1].data
            else:
                imageuncropped = hdulist[0].data
            image = imageuncropped[55:165,60:180]
            
            # display map on zscale
            norm = ImageNormalize(data=image, interval=ZScaleInterval(), stretch=LinearStretch())
            im = ax[a,b].pcolorfast(image, norm=norm, cmap='viridis')
            kernel = Gaussian2DKernel(x_stddev=9)
            convimage = convolve(image, kernel)
            ax[a,b].contour(convimage, levels=np.linspace(convimage.min(), convimage.max(), 12), colors='whitesmoke', alpha=0.7, origin='lower', linewidths=1, zorder=1)
            
            # Species annotation
            props = dict(boxstyle='square', facecolor='k', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(8, 100, 'CN', color='white',fontsize=12, fontweight='bold',bbox=props)
            
            # Date annotation
            ax[a,b].text(80, 100, date, color='white',fontsize=12, fontweight='bold',bbox=props)
            
            # Get average angle values for the date
            table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
            names = table['new_arcfile_name']
            rows = table[names.str.contains(date)]
            velo = rows['velocityPA']
            sun = rows['sunTargetPA']
            avvelo = 360-np.mean(velo) #degrees clockwise from north
            avsun = 360-np.mean(sun) #degrees clockwise from north
            
            # antisunward arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avsun* np.pi / 180)
            dy = l * np.cos(avsun* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='dimgrey', ec='k',head_width=7,head_length=7,width=2, length_includes_head=True)
            # -S annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            props = dict(boxstyle='square', facecolor='black', alpha=0.7, pad=0.2, edgecolor='k')
            ax[a,b].text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=12, fontweight='bold', bbox=props) #annotation for arrow
            
            # antivelocity arrow
            start = (5,5)
            l = 30
            dx = l * np.sin(avvelo* np.pi / 180)
            dy = l * np.cos(avvelo* np.pi / 180) 
            ax[a,b].arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=7,head_length=7,width=2,length_includes_head=True)
            # -v annotation
            txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
            ax[a,b].text(txtpos[0], txtpos[1], '-V', color='white',fontsize=10, fontweight='bold', bbox = props) #annotation for arrow
        
    
            #Scale annotation
            dates = ['2019-11-14', '2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
            scale_date = [31.0997061, 32.9682025, 34.209806268718864, 34.328958743007014, 35.47341127722916, 35.532706858016816, 35.59324266213401, 35.57226744800711, 33.07054419420382, 32.84287518811179, 31.377380152136006, 30.253009538172453, 29.87299630886403, 27.75608220392316, 27.38499010989345] #from table of 10kkm apertures
            width = scale_date[dates.index(date)]
            centrepoint = 95
            ax[a,b].arrow(centrepoint-width/2, 10, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=2,length_includes_head=True)
            ax[a,b].text(81, 14, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
            
            ax[a,b].axis('off')
            ax[a,b].set_aspect('equal')
     
    ax[0, 0].axis('off')
    ax[0, 1].axis('off')
    ax[1, 0].axis('off')
    ax[1, 1].axis('off')
    plt.tight_layout(pad=0, w_pad=0, h_pad=0)
    ax[0, 0].remove()
    ax[0, 1].remove()
    ax[1, 0].remove()
    ax[1, 1].remove()
    
    bbox = fig.bbox_inches.from_bounds(-0.2, 0, 19,10.5)
    plt.savefig('../output_images/Compare_allflux_CN_elongation_nov.pdf', bbox_inches=bbox)
    
    return
            
            
grid_allflux_nov()  
grid_allflux_dec()  
grid_allflux_feb()      
grid_allflux_march()     