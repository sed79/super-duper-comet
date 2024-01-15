#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:10:37 2023

@author: sed79

Gets median S' value under 5,000km diameter aperture over 2I from the 2D arrays of 
solar reflectance values. Finds root mean squared error for uncertainties with

S'uncertainty = sqrt( (S'_i - S'median)**2 / N ) 

"""

import numpy as np

dates = ['2019-11-14', '2019-11-15', '2019-11-26', '2019-12-05', 
         '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', 
         '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', 
         '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
scale_date = [31.0997061, 31.0997061, 32.9682025, 34.209806268718864, 34.328958743007014, 35.47341127722916, 35.532706858016816, 35.59324266213401, 35.57226744800711, 33.07054419420382, 32.84287518811179, 31.377380152136006, 30.253009538172453, 29.87299630886403, 27.75608220392316, 27.38499010989345] #from table of 10kkm apertures
        
def median():
    S_values = []
    S_uncertanties = []
    
    for date in dates:
        flux_vals = []
        num_pix = 0
        sum_square_variances = 0
        grads = np.load('../data/MapsForSlopes/ReflectivityGradients_'+date+'.npy')
        tenthousandkm_in_pixels = scale_date[dates.index(date)]
        for x in range(grads.shape[0]):
            for y in range(grads.shape[1]):
                if np.sqrt((x-200)**2+(y-200)**2) <= (tenthousandkm_in_pixels/4): #aperture diameter of 5,000km
                    flux_vals.append(grads[x, y])
                    num_pix +=1
        median_S = np.median(flux_vals)
        S_values.append(median_S)
        for x in range(grads.shape[0]):
            for y in range(grads.shape[1]):
                if np.sqrt((x-200)**2+(y-200)**2) <= (tenthousandkm_in_pixels/4): #aperture diameter of 5,000km
                    var = grads[x, y] - median_S
                    sum_square_variances += (var**2)
        S_uncertainy = np.sqrt(sum_square_variances / num_pix)
        S_uncertanties.append(S_uncertainy)
    return S_values, S_uncertanties

#final values are saved in extra/MUSES'forPlotting.csv
        

    
