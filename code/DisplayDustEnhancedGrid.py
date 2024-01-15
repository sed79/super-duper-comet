#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:29:24 2023

@author: sed79

Displays a grid of the enhanced dust maps.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from astropy.visualization import simple_norm, ZScaleInterval, LinearStretch, ImageNormalize, ManualInterval, PercentileInterval
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve

        


directory =  '../data/Enhanced_200200/' #200,200 centered maps

dates_undiff = ['2019-11-14', '2019-11-15', '2019-11-26', '2019-12-05', 
         '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29']
dates_diff = ['2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', 
         '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
dates  = ['2019-11-14', '2019-11-15', '2019-11-26', '2019-12-05', 
         '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29',
         '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', 
         '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']

def plot_DUST(to_plot):
    fig, axs = plt.subplots(4,4,figsize=(13,13), gridspec_kw={'height_ratios':[1, 1, 1, 1.15]})
    [(ax1, ax2, ax3, ax4) , (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12), (ax13, ax14, ax15, ax16)] = axs
    axes = [ax1, ax2, ax3, ax4 , ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16]
    
    props = dict(boxstyle='square', facecolor='black', alpha=0.7, pad=0.2, edgecolor='k')
    stddev = 1
    numlevels = 8
    
    for ax, date in zip(axes,dates):
        # diff maps and original stored in different locations
        if date in dates_undiff:
            filename = 'DUST_'+date+'.'+to_plot+'.fits'
        elif date in dates_diff:
            filename = 'DUST_'+date+'_diff.'+to_plot+'.fits'
            
        hdulist = fits.open(directory+filename,ignore_missing_end=True)
        imageuncrop = hdulist[0].data #400,400 pixel image
        # imagecrop = imageuncrop[100:340, 90:330] # this cropping step is previously done in another script
        # image = imagecrop[55:165,60:180]
        image = imageuncrop[115:285, 110:290] # this cropping step is previously done in another script
        kernel = Gaussian2DKernel(x_stddev=stddev)
        convimage = convolve(image, kernel)
    
        #norm = ImageNormalize(data=image, interval=PercentileInterval(99.5), stretch=LinearStretch())
        norm = ImageNormalize(data=image, vmin=-300, vmax=300, stretch=LinearStretch())
        im = ax.pcolorfast(image, norm=norm, cmap='gist_heat')
        
        # Center cross
        #ax.scatter(90,85, marker='x', color='k')

        # Date annotation
        ax.text(110,155, date, color='white',fontsize=12, fontweight='bold',bbox=props)
        
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
        l = 45
        dx = l * np.sin(avsun* np.pi / 180)
        dy = l * np.cos(avsun* np.pi / 180) 
        ax.arrow(start[0], start[1], dx, dy,fc='dimgrey', ec='k',head_width=9,head_length=9,width=3,length_includes_head=True)
        # -S annotation
        txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
        ax.text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=12, fontweight='bold', bbox=props) #annotation for arrow
        
        # antivelocity arrow
        start = (5,5)
        l = 45
        dx = l * np.sin(avvelo* np.pi / 180)
        dy = l * np.cos(avvelo* np.pi / 180) 
        ax.arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=9,head_length=9,width=3,length_includes_head=True)
        # -v annotation
        txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
        ax.text(txtpos[0], txtpos[1], '-V', color='white',fontsize=10, fontweight='bold', bbox = props) #annotation for arrow
    

        #Scale annotation
        #dates = ['2019-11-14', '2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
        scale_date = [31.0997061, 31.0997061, 32.9682025, 34.209806268718864, 34.328958743007014, 35.47341127722916, 35.532706858016816, 35.59324266213401, 35.57226744800711, 33.07054419420382, 32.84287518811179, 31.377380152136006, 30.253009538172453, 29.87299630886403, 27.75608220392316, 27.38499010989345] #from table of 10kkm apertures
        width = scale_date[dates.index(date)]
        centrepoint = 145
        ax.arrow(centrepoint-width/2, 12, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=3,length_includes_head=True)
        ax.text(120, 20, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
        
        ax.axis('off')
        ax.set_aspect('equal')
        
        plt.tight_layout()

    fig.subplots_adjust(bottom=0.05, wspace=0.05, hspace=0.05)
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='gist_heat'),ax=(ax13, ax14, ax15, ax16), extend='both', label="Enhancement Value", fraction=0.1, shrink=0.8, location='bottom', pad=0.04)
    
        
    plt.savefig('../output_images/EnhancedDustGrid_diff.pdf')
        
    return 
        
to_plot = 'azmedsub'
#to_plot = 'azmeddiv'
plot_DUST(to_plot)        
