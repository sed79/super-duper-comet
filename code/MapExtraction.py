#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 11:20:20 2023

@author: sed79

This will create C2, and NH2 maps from dust subtracted cubes, CN maps from 
CN subtracted cubes, and dust maps from the shifted cubes.
"""

import os
import numpy as np
from astropy.io import fits
from mpdaf.obj import Cube
import astropy.units as u

ranges1 = [(5075, 5171)] #C2
ranges2 = [(5973,5982), (5991,6001), (6016,6026), (6330,6340)] #NH2
ranges3 = [(9135,9210)] #CN
ranges4 = [(7080,7120)] #dust

def make_C2_NH2():
    directory = '../data/DustSubtracted/'
    for i, filename in enumerate(sorted(os.listdir(directory))):
        cube = Cube(directory+filename)
        
        for region in ranges1:
            subcube = cube.select_lambda(lbda_min=region[0], lbda_max=region[1], unit_wave=u.angstrom)
            image = subcube.sum(axis=0)
            writing_name = 'C2_'+filename[15:]
            hdu = fits.PrimaryHDU(image.data.data)
            hdu.writeto('../data/ExtractedMaps/C2/'+writing_name)

        NH2 = []
        for region in ranges2:
            subcube = cube.select_lambda(lbda_min=region[0], lbda_max=region[1], unit_wave=u.angstrom)
            image = subcube.sum(axis=0)
            NH2.append(image)
            
        image1, image2, image3, image4 = NH2
        image = image1+image2+image3+image4

        writing_name = 'NH2_'+filename[15:] 
        hdu = fits.PrimaryHDU(image.data.data)
        hdu.writeto('../data/ExtractedMaps/NH2/'+writing_name)
    return

def make_CN():
    directory = '../data/CNSubtracted/'
    for i, filename in enumerate(sorted(os.listdir(directory))):
        cube = Cube(directory+filename)
        
        for region in ranges3:
            subcube = cube.select_lambda(lbda_min=region[0], lbda_max=region[1], unit_wave=u.angstrom)
            image = subcube.sum(axis=0)
            writing_name = 'CN_'+filename[6:]
            hdu = fits.PrimaryHDU(image.data.data)
            hdu.writeto('../data/ExtractedMaps/CN/'+writing_name)
    return
            
def make_DUST():
    directory = '../data/Shifted/'
    for i, filename in enumerate(sorted(os.listdir(directory))):
        cube = Cube(directory+filename)
        
        for region in ranges4:
            subcube = cube.select_lambda(lbda_min=region[0], lbda_max=region[1], unit_wave=u.angstrom)
            image = subcube.sum(axis=0)
            writing_name = 'DUST_'+filename[8:]
            hdu = fits.PrimaryHDU(image.data.data)
            hdu.writeto('../data/ExtractedMaps/DUST/'+writing_name)    
    return
           