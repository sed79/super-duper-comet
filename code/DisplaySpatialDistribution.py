#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 15:35:37 2023

@author: sed79
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.text as mtext
from astropy.io import fits
from astropy.visualization import simple_norm, ZScaleInterval, LinearStretch, ImageNormalize, ManualInterval, PercentileInterval
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.visualization import make_lupton_rgb
from matplotlib.patches import Rectangle

def option1():
    """
    Individual maps zmax
    Colour map zmax

    Returns
    -------
    None.

    """
    directory = '../data/CroppedMaps/EarlyDecTrio/'
    date = '2019-12-'
    fig, ax = plt.subplots(1, 4, figsize=(14, 4.5), gridspec_kw={'width_ratios': [0.95, 1, 1, 1]})
    images = []
    species = ['NH$_2$','C$_2$',  'CN']
    colours = ['red', 'lime', 'xkcd:vivid blue']
    manual_locations = [[(110, 19), (110, 50), (110,73), (110, 90)],
                        [(110, 0), (110, 20), (110,50), (115, 78)],
                        [(110, 20), (110, 47), (110,75), (110, 90)]]
    labels = [['20%', '40%', '60%', '80%'],
              ['', '40%', '60%', '80%'],
              ['20%', '40%', '60%', '80%']]
    props_clabel = dict(boxstyle='round4', facecolor='k', alpha=0.5, pad=0.2, edgecolor='None')
    for i, filename in enumerate(sorted(os.listdir(directory))):
            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            image = hdulist[0].data
            numlevels = 6
            stddev = 5
            kernel = Gaussian2DKernel(x_stddev=stddev)
            convimage = convolve(image, kernel)
            
            #histogram
            ax_histx = ax[i+1].inset_axes([0, 1.05, 1, 0.25], sharex=ax[i+1])
            ax_histx.tick_params(axis="x", labelbottom=False, bottom=False)
            #x_flux = np.sum(image, axis=0) 
            x_flux = image[100]
            ax_histx.bar(range(image.shape[1]), x_flux, width=1, color='dimgrey')
            ax_histx.set(ylim=(0, x_flux.max()))
            
            # Get average angle values for the month december
            table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
            names = table['new_arcfile_name']
            rows = table[names.str.contains(date)]
            velo = rows['velocityPA']
            sun = rows['sunTargetPA']
            avvelo = 360-np.mean(velo) #degrees clockwise from north
            avsun = 360-np.mean(sun) #degrees clockwise from north
            
            # Species annotation
            props = dict(boxstyle='square', facecolor='k', alpha=0.5, pad=0.2, edgecolor='None')
            ax[i+1].text(20, 210, species[i], color='white',fontsize=14, fontweight='bold',bbox=props)
            
            # display individual map
            norm = ImageNormalize(data=image, interval=ZScaleInterval(), vmax=image.max(), stretch=LinearStretch()) #
            im = ax[i+1].pcolor(image, norm=norm, cmap='Greys_r')
            cbar = fig.colorbar(im, ax=ax[i+1], fraction=0.1, shrink=0.9, location='bottom', pad=0.02)
            
            levels = np.linspace(convimage.min(), convimage.max(), numlevels, endpoint=True)
            #print(levels)
            #print(((levels-levels.min()) / (levels.max()-levels.min()) * 100)) #check the flux levels are 20% gaps of conv flux
            #print(((levels) / (image.max()) * 100))
            cont = ax[i+1].contour(convimage, levels=levels, colors=colours[i], origin='lower', zorder=1)
            #cbar.add_lines(cont)
            tl = ax[i+1].clabel(cont, cont.levels, inline=False, inline_spacing=0.5, fontsize=7, colors='w',use_clabeltext=True, manual=manual_locations[i])
            for j, t in enumerate(tl):
                t.set_bbox(props_clabel)
                t.set_text(labels[i][j])
                t.set_rotation('horizontal')
                t.set_fontweight('bold')
                
            

            #set coloured image scale
            images.append(norm(image))
            
    # make and display final rgb image       
    image_r, image_g, image_b = images #c2, nh2, cn
    rgbimage = make_lupton_rgb(image_r, image_g, image_b, stretch=1, Q=0)   
    im = ax[0].pcolorfast(rgbimage)
    ax[0].text(20, 210, 'r: NH$_2$', color='white',fontsize=14, fontweight='bold',bbox=props)
    ax[0].text(20, 185, 'g: C$_2$', color='white',fontsize=14, fontweight='bold',bbox=props)
    ax[0].text(20, 160, 'b: CN', color='white',fontsize=14, fontweight='bold',bbox=props)
    
    
    # antisunward arrow
    start = (10, 10)
    l = 50
    dx = l * np.sin(avsun* np.pi / 180)
    dy = l * np.cos(avsun* np.pi / 180) 
    ax[0].arrow(start[0], start[1], dx, dy,fc='xkcd:fire engine red', ec='k',head_width=10,head_length=10,width=4,length_includes_head=True)
    # -S annotation
    txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
    props = dict(boxstyle='square', facecolor='black', alpha=0.5, pad=0.2, edgecolor='None')
    ax[0].text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=12, fontweight='bold', bbox=props) #annotation for arrow
    # antivelocity arrow
    start = (10, 10)
    l = 50
    dx = l * np.sin(avvelo* np.pi / 180)
    dy = l * np.cos(avvelo* np.pi / 180) 
    ax[0].arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=10,head_length=10,width=4,length_includes_head=True)
    # -v annotation
    txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
    ax[0].text(txtpos[0], txtpos[1], '-V', color='white',fontsize=10, fontweight='bold', bbox = props) #annotation for arrow
    # Coordinate arrows
    start = (230, 10)
    ax[0].arrow(start[0], start[1], 0, 50,fc='k', ec='k',head_width=10,head_length=10,width=3,length_includes_head=True)
    ax[0].arrow(start[0]+1, start[1], -50, 0,fc='k', ec='k',head_width=10,head_length=10,width=3,length_includes_head=True)
    # Coordinate text
    ax[0].text(225, 68, 'N', color='white',fontsize=10, fontweight='bold', bbox=props) #inline with arrows with small adjustments
    ax[0].text(167, 6, 'E', color='white',fontsize=10, fontweight='bold', bbox=props) 
    
    #Scale annotation
    scalesdec = [34.209806268718864, 34.328958743007014]
    width = np.mean(scalesdec)
    centrepoint = 180
    ax[0].arrow(centrepoint-width/2, 200, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=4,length_includes_head=True)
    ax[0].text(150, 210, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
    for i in range(4):
        ax[i].axis('off')
        ax[i].set_aspect('equal')
    
    
    plt.tight_layout()
    #bbox = fig.bbox_inches.from_bounds(0, 1, 14, 4)
    plt.savefig('../output_images/Colour_trio_contours_hist_option1.pdf')#, bbox_inches=bbox)

    return

#option1()

def option2():
    """
    Individual maps zmax
    Colour map zmax

    Returns
    -------
    None.

    """
    directory = '../data/CroppedMaps/EarlyDecTrio/'
    date = '2019-12-'
    fig, ax = plt.subplots(4, 1, figsize=(6, 14))
    images = []
    species = ['NH$_2$','C$_2$',  'CN']
    colours = ['white', 'white', 'white']
    manual_locations = [[(110, 19), (110, 50), (110,73), (110, 90)],
                        [(110, 0), (110, 20), (110,50), (115, 78)],
                        [(110, 20), (110, 47), (110,75), (110, 90)]]
    labels = [['20%', '40%', '60%', '80%'],
              ['', '40%', '60%', '80%'],
              ['20%', '40%', '60%', '80%']]
    props_clabel = dict(boxstyle='round4', facecolor='k', alpha=0.5, pad=0.2, edgecolor='None')
    props = dict(boxstyle='square', facecolor='k', alpha=0.7, pad=0.2, edgecolor='None')
    for i, filename in enumerate(sorted(os.listdir(directory))):
            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            image = hdulist[0].data
            numlevels = 6
            stddev = 5
            kernel = Gaussian2DKernel(x_stddev=stddev)
            convimage = convolve(image, kernel)
            
            #histogram x
            ax_histx = ax[i+1].inset_axes([0, 1, 1, 0.15], sharex=ax[i+1]) #[posx, posy, width, heigth]
            ax_histx.tick_params(axis="x", labelbottom=False, bottom=False)
            #x_flux = np.sum(image, axis=0) #sum all flux up image to display in histogram
            horizontal_strip = image[95:105] #select area of flux +- 5 pixels around comet
            x_flux = np.sum(horizontal_strip, axis=0) #y=100
            ax_histx.bar(range(image.shape[1]), x_flux, width=1, color='dimgrey')
            ax_histx.set(ylim=(0, x_flux.max()), xlim=(0,240))
            
            
            #histogram y
            ax_histy = ax[i+1].inset_axes([1, 0, 0.15, 1], sharey=ax[i+1])
            ax_histy.tick_params(axis="y", labelleft=False, left=False)
            #y_flux = np.sum(image, axis=1) 
            vertical_strip = image[:,105:115] #select area of flux +- 5 pixels around comet
            y_flux = np.sum(vertical_strip, axis=1) #y=100
            #y_flux = image[:,110] #x=110
            ax_histy.barh(np.arange(image.shape[0]), y_flux, color='dimgrey')
            ax_histy.set(xlim=(0, y_flux.max()), ylim=(0,240))
            
            # Dashed lines
            ax[i+1].vlines(105, ymin=0, ymax=240, linestyle='dashed', colors='white')
            ax[i+1].vlines(115, ymin=0, ymax=240, linestyle='dashed', colors='white')
            ax[i+1].hlines(100, xmin=0, xmax=240, linestyle='dashed', colors='white')
            
            # Get average angle values for the month december
            table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
            names = table['new_arcfile_name']
            rows = table[names.str.contains(date)]
            velo = rows['velocityPA']
            sun = rows['sunTargetPA']
            avvelo = 360-np.mean(velo) #degrees clockwise from north
            avsun = 360-np.mean(sun) #degrees clockwise from north
            
            # Species annotation
            ax[i+1].text(20, 210, species[i], color='white',fontsize=14, fontweight='bold',bbox=props)
            
            # display individual map
            norm = ImageNormalize(data=image, interval=ZScaleInterval(), vmax=image.max(), stretch=LinearStretch()) #
            im = ax[i+1].pcolorfast(image, norm=norm, cmap='viridis')
            cbar = fig.colorbar(im, ax=ax[i+1], fraction=0.1, shrink=0.9, location='left', pad=0.02)
            
            levels = np.linspace(convimage.min(), convimage.max(), numlevels, endpoint=True)
            #print(levels)
            #print(((levels-levels.min()) / (levels.max()-levels.min()) * 100)) #check the flux levels are 20% gaps of conv flux
            #print(((levels) / (image.max()) * 100))
            cont = ax[i+1].contour(convimage, levels=levels, colors=colours[i], origin='lower', zorder=1)
            #cbar.add_lines(cont)
            tl = ax[i+1].clabel(cont, cont.levels, inline=False, inline_spacing=0.5, fontsize=7, colors='w',use_clabeltext=True, manual=manual_locations[i])
            for j, t in enumerate(tl):
                t.set_bbox(props_clabel)
                t.set_text(labels[i][j])
                t.set_rotation('horizontal')
                t.set_fontweight('bold')
                
            

            #set coloured image scale
            images.append(norm(image))
            
    # make and display final rgb image       
    image_r, image_g, image_b = images #c2, nh2, cn
    rgbimage = make_lupton_rgb(image_r, image_g, image_b, stretch=1, Q=0)   
    im = ax[0].pcolorfast(rgbimage)
    ax[0].text(20, 210, 'r: NH$_2$', color='white',fontsize=14, fontweight='bold',bbox=props)
    ax[0].text(20, 185, 'g: C$_2$', color='white',fontsize=14, fontweight='bold',bbox=props)
    ax[0].text(20, 160, 'b: CN', color='white',fontsize=14, fontweight='bold',bbox=props)
    
    
    # antisunward arrow
    start = (10, 10)
    l = 50
    dx = l * np.sin(avsun* np.pi / 180)
    dy = l * np.cos(avsun* np.pi / 180) 
    ax[0].arrow(start[0], start[1], dx, dy,fc='xkcd:fire engine red', ec='k',head_width=10,head_length=10,width=4,length_includes_head=True)
    # -S annotation
    txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
    ax[0].text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=12, fontweight='bold', bbox=props) #annotation for arrow
    # antivelocity arrow
    start = (10, 10)
    l = 50
    dx = l * np.sin(avvelo* np.pi / 180)
    dy = l * np.cos(avvelo* np.pi / 180) 
    ax[0].arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=10,head_length=10,width=4,length_includes_head=True)
    # -v annotation
    txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
    ax[0].text(txtpos[0], txtpos[1], '-V', color='white',fontsize=10, fontweight='bold', bbox = props) #annotation for arrow
    # Coordinate arrows
    start = (230, 10)
    ax[0].arrow(start[0], start[1], 0, 50,fc='k', ec='k',head_width=10,head_length=10,width=3,length_includes_head=True)
    ax[0].arrow(start[0]+1, start[1], -50, 0,fc='k', ec='k',head_width=10,head_length=10,width=3,length_includes_head=True)
    # Coordinate text
    ax[0].text(225, 68, 'N', color='white',fontsize=10, fontweight='bold', bbox=props) #inline with arrows with small adjustments
    ax[0].text(167, 6, 'E', color='white',fontsize=10, fontweight='bold', bbox=props) 
    
    #Scale annotation
    scalesdec = [34.209806268718864, 34.328958743007014]
    width = np.mean(scalesdec)
    centrepoint = 180
    ax[0].arrow(centrepoint-width/2, 200, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=4,length_includes_head=True)
    ax[0].text(150, 210, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
    for i in range(4):
        ax[i].axis('off')
        ax[i].set_aspect('equal')
        
    
    
    plt.tight_layout()
    #bbox = fig.bbox_inches.from_bounds(0, 1, 14, 4)
    plt.savefig('../output_images/Colour_trio_contours_hist_option2.pdf')#, bbox_inches=bbox)

    return

def option3():
    """
    Individual maps zmax
    Colour map zmax

    Returns
    -------
    None.

    """
    directory = '../data/CroppedMaps/EarlyDecTrio/'
    date = '2019-12-'
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    ax1, ax2, ax3 = axes
    axs = [ax2, ax1, ax3] #rearranging images without changing file structure
    images = []
    species = ['NH$_2$','C$_2$',  'CN']
    colours = ['white', 'white', 'white']
    manual_locations = [[(190, 25), (150, 60), (130,80), (125, 90)],
                        [(110, 0), (200, 30), (145,55), (125, 85)],
                        [(215, 25), (150, 50), (130,70), (120, 93)]]
    labels = [['20%', '40%', '60%', '80%'],
              ['', '40%', '60%', '80%'],
              ['20%', '40%', '60%', '80%']]
    props_clabel = dict(boxstyle='round4', facecolor='k', alpha=0.7, pad=0.2, edgecolor='None')
    props = dict(boxstyle='square', facecolor='k', alpha=0.9, pad=0.2, edgecolor='None')
    for i, (ax, filename) in enumerate(zip(axs,sorted(os.listdir(directory)))):
            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            image = hdulist[0].data
            numlevels = 6
            stddev = 5
            kernel = Gaussian2DKernel(x_stddev=stddev)
            convimage = convolve(image, kernel)
            
            #histogram x
            ax_histx = ax.inset_axes([0, 1, 1, 0.15], sharex=ax) #[posx, posy, width, heigth]
            ax_histx.tick_params(axis="x", labelbottom=False, bottom=False)
            #x_flux = np.sum(image, axis=0) #sum all flux up image to display in histogram
            horizontal_strip = image[95:105] #select area of flux +- 5 pixels around comet
            x_flux = np.sum(horizontal_strip, axis=0) #y=100
            ax_histx.bar(range(image.shape[1]), x_flux, width=1, color='dimgrey')
            ax_histx.set(ylim=(0, x_flux.max()), xlim=(0,240))
            
            
            #histogram y
            ax_histy = ax.inset_axes([1, 0, 0.15, 1], sharey=ax)
            ax_histy.tick_params(axis="y", labelleft=False, left=False)
            #y_flux = np.sum(image, axis=1) 
            vertical_strip = image[:,105:115] #select area of flux +- 5 pixels around comet
            y_flux = np.sum(vertical_strip, axis=1) #y=100
            #y_flux = image[:,110] #x=110
            ax_histy.barh(np.arange(image.shape[0]), y_flux, color='dimgrey')
            ax_histy.set(xlim=(0, y_flux.max()), ylim=(0,240))
            
            # Dashed lines
            if i == 0:
                ax.add_patch(Rectangle((0,95), 240, 10, alpha=0.5,facecolor='grey', edgecolor='None')) #horizontal
                ax.add_patch(Rectangle((105,105), 10, 135, alpha=0.5,facecolor='grey', edgecolor='None')) #vertical top half
                ax.add_patch(Rectangle((105,0), 10, 95, alpha=0.5,facecolor='grey', edgecolor='None')) #vertical bottom half
            elif i == 1:
                ax.vlines(105, ymin=0, ymax=95, linestyle='dashed', colors='white') #lower
                ax.vlines(115, ymin=0, ymax=95, linestyle='dashed', colors='white') #lower
                ax.vlines(105, ymin=105, ymax=240, linestyle='dashed', colors='white') #upper
                ax.vlines(115, ymin=105, ymax=240, linestyle='dashed', colors='white')
                ax.hlines(95, xmin=0, xmax=105, linestyle='dashed', colors='white') #left
                ax.hlines(105, xmin=0, xmax=105, linestyle='dashed', colors='white')
                ax.hlines(95, xmin=115, xmax=240, linestyle='dashed', colors='white') #right
                ax.hlines(105, xmin=115, xmax=240, linestyle='dashed', colors='white')
            elif i == 2:
                ax.vlines(105, ymin=0, ymax=95, linestyle='dashed', colors='grey') #lower
                ax.vlines(115, ymin=0, ymax=95, linestyle='dashed', colors='grey') #lower
                ax.vlines(105, ymin=105, ymax=240, linestyle='dashed', colors='grey') #upper
                ax.vlines(115, ymin=105, ymax=240, linestyle='dashed', colors='grey')
                ax.hlines(95, xmin=0, xmax=105, linestyle='dashed', colors='grey') #left
                ax.hlines(105, xmin=0, xmax=105, linestyle='dashed', colors='grey')
                ax.hlines(95, xmin=115, xmax=240, linestyle='dashed', colors='grey') #right
                ax.hlines(105, xmin=115, xmax=240, linestyle='dashed', colors='grey')
            
            # Get average angle values for the month december
            table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
            names = table['new_arcfile_name']
            rows = table[names.str.contains(date)]
            velo = rows['velocityPA']
            sun = rows['sunTargetPA']
            avvelo = 360-np.mean(velo) #degrees clockwise from north
            avsun = 360-np.mean(sun) #degrees clockwise from north
            
            # Species annotation
            ax.text(20, 210, species[i], color='white',fontsize=14, fontweight='bold',bbox=props)
            
            # display individual map
            norm = ImageNormalize(data=image, interval=ZScaleInterval(), vmax=image.max(), stretch=LinearStretch()) #
            im = ax.pcolorfast(image, norm=norm, cmap='viridis')
            cbar = fig.colorbar(im, ax=ax, fraction=0.1, shrink=0.9, location='bottom', pad=0.02)
            
            levels = np.linspace(convimage.min(), convimage.max(), numlevels, endpoint=True)
            #print(levels)
            #print(((levels-levels.min()) / (levels.max()-levels.min()) * 100)) #check the flux levels are 20% gaps of conv flux
            #print(((levels) / (image.max()) * 100))
            cont = ax.contour(convimage, levels=levels, colors=colours[i], origin='lower', zorder=1)
            #cbar.add_lines(cont)
            tl = ax.clabel(cont, cont.levels, inline=False, inline_spacing=0.5, fontsize=7, colors='w',use_clabeltext=True, manual=manual_locations[i])
            for j, t in enumerate(tl):
                t.set_bbox(props_clabel)
                t.set_text(labels[i][j])
                t.set_rotation('horizontal')
                t.set_fontweight('bold')
                
            ax.axis('off')
            ax.set_aspect('equal')

            #set coloured image scale
            images.append(norm(image))
            
    # antisunward arrow
    start = (10, 10)
    l = 50
    dx = l * np.sin(avsun* np.pi / 180)
    dy = l * np.cos(avsun* np.pi / 180) 
    ax1.arrow(start[0], start[1], dx, dy,fc='dimgrey', ec='k',head_width=10,head_length=10,width=4,length_includes_head=True)
    # -S annotation
    txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
    ax1.text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=12, fontweight='bold', bbox=props) #annotation for arrow
    # antivelocity arrow
    start = (10, 10)
    l = 50
    dx = l * np.sin(avvelo* np.pi / 180)
    dy = l * np.cos(avvelo* np.pi / 180) 
    ax1.arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=10,head_length=10,width=4,length_includes_head=True)
    # -v annotation
    txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
    ax1.text(txtpos[0], txtpos[1], '-V', color='white',fontsize=10, fontweight='bold', bbox = props) #annotation for arrow
    
    #Scale annotation
    scalesdec = [34.209806268718864, 34.328958743007014]
    width = np.mean(scalesdec)
    centrepoint = 180
    ax1.arrow(centrepoint-width/2, 200, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=3,length_includes_head=True)
    ax1.text(150, 210, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
    
    plt.tight_layout()
    #bbox = fig.bbox_inches.from_bounds(0, 1, 14, 4)
    plt.savefig('../output_images/Colour_trio_contours_hist_option3.pdf')#, bbox_inches=bbox)

    return

def option4():
    """
    Individual maps zmax
    Colour map zmax

    Returns
    -------
    None.

    """
    directory = '../data/CroppedMaps/EarlyDecTrio/'
    date = '2019-12-'
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    ax1, ax2, ax3 = axes
    axs = [ax2, ax1, ax3] #rearranging images without changing file structure
    images = []
    species = ['NH$_2$','C$_2$',  'CN']
    colours = ['white', 'white', 'white']
    manual_locations = [[(110, 19), (110, 50), (110,73), (110, 90)],
                        [(110, 0), (110, 20), (110,40), (115, 78)],

                        [(110, 20), (110, 47), (110,75), (110, 90)]]
    labels = [['20%', '40%', '60%', '80%'],
              ['', '40%', '60%', '80%'],
              ['20%', '40%', '60%', '80%']]
    props_clabel = dict(boxstyle='round4', facecolor='k', alpha=0.7, pad=0.2, edgecolor='None')
    props = dict(boxstyle='square', facecolor='k', alpha=0.9, pad=0.2, edgecolor='None')
    for i, (ax, filename) in enumerate(zip(axs,sorted(os.listdir(directory)))):
            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            image = hdulist[0].data
            numlevels = 6
            stddev = 5
            kernel = Gaussian2DKernel(x_stddev=stddev)
            convimage = convolve(image, kernel)
            
            #histogram x
            ax_histx = ax.inset_axes([0, 1, 1, 0.15], sharex=ax) #[posx, posy, width, heigth]
            ax_histx.tick_params(axis="both", labelbottom=False, bottom=False, labelleft=False, left=False)
            #x_flux = np.sum(image, axis=0) #sum all flux up image to display in histogram
            horizontal_strip = image[95:105] #select area of flux +- 5 pixels around comet
            x_flux = np.sum(horizontal_strip, axis=0) #y=100
            ax_histx.bar(range(image.shape[1]), x_flux, width=1, color='dimgrey')
            ax_histx.set(ylim=(0, x_flux.max()), xlim=(0,240))
            
            
            #histogram y
            ax_histy = ax.inset_axes([1, 0, 0.15, 1], sharey=ax)
            ax_histy.tick_params(axis="both", labelleft=False, left=False, labelbottom=False, bottom=False)
            #y_flux = np.sum(image, axis=1) 
            vertical_strip = image[:,105:115] #select area of flux +- 5 pixels around comet
            y_flux = np.sum(vertical_strip, axis=1) #y=100
            #y_flux = image[:,110] #x=110
            ax_histy.barh(np.arange(image.shape[0]), y_flux, color='dimgrey')
            ax_histy.set(xlim=(0, y_flux.max()), ylim=(0,240))
            

            # Get average angle values for the month december
            table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
            names = table['new_arcfile_name']
            rows = table[names.str.contains(date)]
            velo = rows['velocityPA']
            sun = rows['sunTargetPA']
            avvelo = 360-np.mean(velo) #degrees clockwise from north
            avsun = 360-np.mean(sun) #degrees clockwise from north
            
            # Species annotation
            ax.text(20, 210, species[i], color='white',fontsize=14, fontweight='bold',bbox=props)
            
            # display individual map
            norm = ImageNormalize(data=image, interval=ZScaleInterval(), vmax=image.max(), stretch=LinearStretch()) #
            im = ax.pcolorfast(image, norm=norm, cmap='Blues_r')
            cbar = fig.colorbar(im, ax=ax, fraction=0.05, shrink=0.9, location='bottom', pad=0.02)
            
            levels = np.linspace(convimage.min(), convimage.max(), numlevels, endpoint=True)
            #print(levels)
            #print(((levels-levels.min()) / (levels.max()-levels.min()) * 100)) #check the flux levels are 20% gaps of conv flux
            #print(((levels) / (image.max()) * 100))
            cont = ax.contour(convimage, levels=levels, colors='whitesmoke', alpha = 0.7, origin='lower', zorder=1)
            #cbar.add_lines(cont)
            tl = ax.clabel(cont, cont.levels, inline=False, inline_spacing=0.5, fontsize=7, colors='w',use_clabeltext=True, manual=manual_locations[i])
            for j, t in enumerate(tl):
                t.set_bbox(props_clabel)
                t.set_text(labels[i][j])
                t.set_rotation('horizontal')
                t.set_fontweight('bold')
                
            ax.axis('off')
            ax.set_aspect('equal')

            #set coloured image scale
            images.append(norm(image))
            
    # antisunward arrow
    start = (10, 10)
    l = 50
    dx = l * np.sin(avsun* np.pi / 180)
    dy = l * np.cos(avsun* np.pi / 180) 
    ax1.arrow(start[0], start[1], dx, dy,fc='dimgrey', ec='k',head_width=10,head_length=10,width=4,length_includes_head=True)
    # -S annotation
    txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
    ax1.text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=12, fontweight='bold', bbox=props) #annotation for arrow
    # antivelocity arrow
    start = (10, 10)
    l = 50
    dx = l * np.sin(avvelo* np.pi / 180)
    dy = l * np.cos(avvelo* np.pi / 180) 
    ax1.arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=10,head_length=10,width=4,length_includes_head=True)
    # -v annotation
    txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
    ax1.text(txtpos[0], txtpos[1], '-V', color='white',fontsize=10, fontweight='bold', bbox = props) #annotation for arrow
    
    #Scale annotation
    scalesdec = [34.209806268718864, 34.328958743007014]
    width = np.mean(scalesdec)
    centrepoint = 180
    ax1.arrow(centrepoint-width/2, 203, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=3,length_includes_head=True)
    ax1.text(150, 210, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
    
    #plt.tight_layout()
    bbox = fig.bbox_inches.from_bounds(1, 0, 13, 5)
    plt.savefig('../output_images/Colour_trio_contours_hist_blues.pdf', bbox_inches=bbox)

    return

def just_rgb():
    directory = '../data/CroppedMaps/EarlyDecTrio/' #NH2 then C2 then CN
    date = '2019-12-'
    fig, ax = plt.subplots(figsize=(4,4))
    images = []
    props = dict(boxstyle='square', facecolor='k', alpha=0.5, pad=0.2, edgecolor='None')
    for i, filename in enumerate(sorted(os.listdir(directory))):
            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            image = hdulist[0].data

            # Get average angle values for the month december
            table = pd.read_csv(r'../extra/JPL_query_table.csv', delimiter = ',')    
            names = table['new_arcfile_name']
            rows = table[names.str.contains(date)]
            velo = rows['velocityPA']
            sun = rows['sunTargetPA']
            avvelo = 360-np.mean(velo) #degrees clockwise from north
            avsun = 360-np.mean(sun) #degrees clockwise from north
            
            
            # get zscale norm
            norm = ImageNormalize(data=image, interval=ZScaleInterval(), vmax=image.max(), stretch=LinearStretch()) #

            #set coloured image scale
            images.append(norm(image))
            
    # make and display final rgb image       
    image_r, image_g, image_b = images #NH2 then C2 then CN
    rgbimage = make_lupton_rgb(image_r, image_g, image_b, stretch=1, Q=0) #linear strech
    im = ax.pcolorfast(rgbimage)
    ax.text(20, 210, 'r: NH$_2$', color='white',fontsize=10, fontweight='bold',bbox=props)
    ax.text(20, 185, 'g: C$_2$', color='white',fontsize=10, fontweight='bold',bbox=props)
    ax.text(20, 160, 'b: CN', color='white',fontsize=10, fontweight='bold',bbox=props)
    
    ax.set(aspect='equal')
    ax.axis('off')
    
    # antisunward arrow
    start = (10, 10)
    l = 50
    dx = l * np.sin(avsun* np.pi / 180)
    dy = l * np.cos(avsun* np.pi / 180) 
    ax.arrow(start[0], start[1], dx, dy,fc='dimgrey', ec='k',head_width=10,head_length=10,width=4,length_includes_head=True)
    # -S annotation
    txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
    ax.text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=12, fontweight='bold', bbox=props) #annotation for arrow
    # antivelocity arrow
    start = (10, 10)
    l = 50
    dx = l * np.sin(avvelo* np.pi / 180)
    dy = l * np.cos(avvelo* np.pi / 180) 
    ax.arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=10,head_length=10,width=4,length_includes_head=True)
    # -v annotation
    txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
    ax.text(txtpos[0], txtpos[1], '-V', color='white',fontsize=10, fontweight='bold', bbox = props) #annotation for arrow
    
    #Scale annotation
    scalesdec = [34.209806268718864, 34.328958743007014]
    width = np.mean(scalesdec)
    centrepoint = 180
    ax.arrow(centrepoint-width/2, 203, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=3,length_includes_head=True)
    ax.text(150, 210, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
    
    plt.tight_layout()
    plt.savefig('../output_images/RGB.pdf')
    return
    

option1()
option2()
option3()
option4()
just_rgb()


