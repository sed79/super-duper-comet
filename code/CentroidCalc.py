#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 16:25:31 2023

@author: sed79

Calculates the optocenter of 2D .fits file using Photutils 2D Gaussian fit. 
Fed in LambdaSum datacubes, 2D arrays made by summing spectral axis of 
datacubes.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from mpdaf.obj import Cube, Image
from astropy.visualization import simple_norm
import astropy.units as u
from photutils.centroids import centroid_sources, centroid_2dg, centroid_1dg

def raw_cubes():

    directory = '../data/LambdaSum/'
    x_init = [] #need to put list of guesses initially
    y_init = []
       
    
    for i, filename in enumerate(os.listdir(directory)):
        hdulist = fits.open(directory+filename,ignore_missing_end=True)
        image = hdulist[0].data
        x, y = centroid_sources(image, x_init[i], y_init[i], box_size=15,
                                centroid_func=centroid_2dg)
        print(x[0])  
        print(y[0])
        
        plt.figure(figsize=(8, 4))
        plt.imshow(image, origin='lower', interpolation='nearest')
        plt.scatter(x, y, marker='+', s=80, color='red')
        #plt.xlim(x-50, x+50)
        #plt.ylim(y-50, y+50)
    return
    
def test_cen():
    """
    Compare with
    (138.18983183, 132.19778292)
    (171.3836756, 153.19635283)
    for the datacubes 2019-12-05T07:55:01.516 and 2020-03-19T05:23:42.010

    Returns
    -------
    None.

    """
    directory = '../data/LambdaSum/'
    testfilenames = ['LambdaSum_2019-12-05T07:55:01.516.fits',
                     'LambdaSum_2020-03-19T05:23:42.010.fits']
    test_x_init = [140, 175]
    test_y_init = [130, 150]
    
    for i, filename in enumerate(testfilenames):
        hdulist = fits.open(directory+filename,ignore_missing_end=True)
        image = hdulist[0].data
        x, y = centroid_sources(image, test_x_init[i], test_y_init[i], box_size=15,
                                centroid_func=centroid_2dg)
        print(filename+' x: '+str(x))  
        print(filename+' y: '+str(y))
        
    return
