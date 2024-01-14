#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 13:26:53 2022

@author: sed79

Remove extreme values from datacube.

"""
import os
import numpy as np
import numpy.ma as ma
from mpdaf.obj import Cube
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm

directory = '../data/DIFFCUBES/'

for i, filename in enumerate(sorted(os.listdir(directory))):
    writename = filename[18:]
    raw = Cube(directory+filename)
    
    allval = raw.data.compressed()
    
    minimum = np.percentile(allval, 0.000005)
    maximum = np.percentile(allval, 99.999995)
    
    condition = np.where(np.logical_or(raw.data < minimum, raw.data > maximum))
    
    new = raw.copy()
    
    new.mask_selection(condition)
    
    new.write('../data/ExtremesMasked/Masked_'+writename)
