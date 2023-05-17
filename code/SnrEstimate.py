#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 13:00:52 2023

@author: sed79

Makes images to help estimate the signal to noise ratio for each of
the nights once maps were coadded
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

directory = '../data/CroppedMaps/NH2_ALL/'


for i, filename in enumerate(sorted(os.listdir(directory))):
    date = filename[12:22]
    hdulist = fits.open(directory+filename,ignore_missing_end=True)
    image = hdulist[0].data
    if date == '2019-11-14':
        image = hdulist[1].data
    
    hrow = image[100]
    vrow = image[:,110]
    xx = np.linspace(0, 240, 240)
    
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(10, 6))
    ax1, ax2 = ax
    
    ax1.bar(xx, vrow)
    ax1.bar(xx, image[:,109])
    ax1.bar(xx, image[:,111])
    ax2.bar(xx, hrow)
    ax2.bar(xx, image[99])
    ax2.bar(xx, image[101])
    
    ax1.text(200,700, date)
    
    plt.savefig('../output_images/SNR_estimate/SNR_estimate_NH2_'+date+'.pdf')
            
            
            
            