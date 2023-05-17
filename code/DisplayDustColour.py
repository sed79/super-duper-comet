#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 10:42:06 2023

@author: sed79
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

dates = ['2019-11-14', '2019-11-15', '2019-11-26', '2019-12-05', 
         '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', 
         '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', 
         '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']

def plot():
    fig, axs = plt.subplots(4,4,figsize=(15,15), gridspec_kw={'height_ratios':[1, 1, 1, 1.14]})
    [(ax1, ax2, ax3, ax4) , (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12), (ax13, ax14, ax15, ax16)] = axs
    axes = [ax1, ax2, ax3, ax4 , ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16]
    
    cmap = mpl.cm.get_cmap("rainbow").copy()
    bounds = [0.001,2.0,4.0,6.0,8.0, 10,12,14, 16]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')
    cmap.set_over('white')
    cmap.set_under('white')
    
    props = dict(boxstyle='square', facecolor='black', alpha=0.7, pad=0.2, edgecolor='k')
    
    for ax, date in zip(axes,dates):
        gradsuncropped = np.load('../data/SlopeMaps/ReflectivityGradients_'+date+'.npy') #400x400 arrays of S' values (the output of SolarReflectanceSlopes.py)
        gradsonecrop = gradsuncropped[100:340, 90:330]
        grads = gradsonecrop[55:165,60:180]
        
        im = ax.pcolorfast(grads, norm=norm, cmap=cmap)
        
        # Date annotation
        ax.text(80,100, date, color='white',fontsize=12, fontweight='bold',bbox=props)
        
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
        ax.arrow(start[0], start[1], dx, dy,fc='dimgrey', ec='k',head_width=7,head_length=7,width=2,length_includes_head=True)
        # -S annotation
        txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
        ax.text(txtpos[0], txtpos[1], '-$\odot$', color='white',fontsize=12, fontweight='bold', bbox=props) #annotation for arrow
        
        # antivelocity arrow
        start = (5,5)
        l = 30
        dx = l * np.sin(avvelo* np.pi / 180)
        dy = l * np.cos(avvelo* np.pi / 180) 
        ax.arrow(start[0], start[1], dx, dy,fc='white', ec='k',head_width=7,head_length=7,width=2,length_includes_head=True)
        # -v annotation
        txtpos = (start[0]+1.1*dx, start[1]+dy*1.1)
        ax.text(txtpos[0], txtpos[1], '-V', color='white',fontsize=10, fontweight='bold', bbox = props) #annotation for arrow
    

        #Scale annotation
        #dates = ['2019-11-14', '2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
        scale_date = [31.0997061, 31.0997061, 32.9682025, 34.209806268718864, 34.328958743007014, 35.47341127722916, 35.532706858016816, 35.59324266213401, 35.57226744800711, 33.07054419420382, 32.84287518811179, 31.377380152136006, 30.253009538172453, 29.87299630886403, 27.75608220392316, 27.38499010989345] #from table of 10kkm apertures
        width = scale_date[dates.index(date)]
        centrepoint = 95
        ax.arrow(centrepoint-width/2, 10, width, 0, fc='w', ec='k', head_width=0,head_length=0,width=2,length_includes_head=True)
        ax.text(81, 14, "10,000 km", color='white',fontsize=10, fontweight='bold', bbox=props)
        
        #ax.axis('off')
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        ax.set_aspect('equal')
        
        plt.tight_layout()    
    
    fig.subplots_adjust(bottom=0.05, wspace=0.05, hspace=0.05)
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),ax=(ax13, ax14, ax15, ax16), label="(% / 1000$\AA$)", fraction=0.1, shrink=0.8, location='bottom', pad=0.04)
    
    plt.savefig('../output_images/DustColourGrid.pdf')
        
    return 


plot()

    