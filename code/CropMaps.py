#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 11:33:58 2023

@author: sed79

Crops maps to shape 240x240 with comet centre at (110.5,100.5) after median coadding so borders aren't involved with display
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import simple_norm

def crop():
    directory = '../data/MedianCoaddedMaps/' # place for median co-added maps
    for i, filename in enumerate(sorted(os.listdir(directory))):
        hdulist = fits.open(directory+filename,ignore_missing_end=True)
        scidata = hdulist[0].data
        cropdata = scidata[100:340, 90:330]
        print('cropping '+filename+' from '+str(np.shape(scidata))+' to '+str(np.shape(cropdata)))
        
        fig, ax = plt.subplots()
        norm = simple_norm(cropdata, 'linear', percent=100)
        ax.imshow(cropdata, norm=norm, origin='lower', cmap='viridis')
        plt.title(filename)
        
        hdu = fits.PrimaryHDU(cropdata)
        writing_name = 'Cropped_'+filename

        hdu.writeto('../data/CroppedMaps/'+writing_name)
        
crop()
