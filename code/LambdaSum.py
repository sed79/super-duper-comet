#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 10:40:43 2022

@author: sed79

Creates 2D arrays and saves them as .fits files from summing 
the spectral axis of a datacube.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from mpdaf.obj import Cube, Image
from astropy.visualization import simple_norm
import astropy.units as u

directory = '../data/CUBES/' #where reduced datacubes are kept
    
for filename in os.listdir(directory):
    cube = Cube(directory+filename)
    
    image = cube.sum(axis=0)
    writing_name = 'LambdaSum_'+filename[18:]
    hdu = fits.PrimaryHDU(image.data.data)
    hdu.writeto('../data/LambdaSum/'+writing_name)
