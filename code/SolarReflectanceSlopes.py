#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 08:49:05 2023

@author: sed79

Using reflectance method outlined in paper about C/2015 ER61 DOI:10.1051/0004-6361/201833960. 

First, a datacube is divided by the solar spectrum which is interpolated to match 
the datacube's spectral coordinates. Then, three wavelenth ranges are extracted 
from the datacube and collapsed along the spectral dimension by taking an average 
to get three maps. The three maps' wavelengths ranges are:

bluemap = 5260A-5300A centered at 5280A
redmap = 7060A-7100A centered at 7080A
infraredmap = 8580A-8620A centered at 8600A

Each map from the same observing epoch are then median coadded so there is one 
blue, red, and infrared map per night.

The normalised reflectivity gradient, S', is then calculated for each pixel using the formula 

S' = ((S2-S1)/(8600-5280)) * (1/meanS) * 100000

where S2 and S1 are the pixel flux values from the infraredmap and blue maps, meanS is 
the flux value from the red map, and 8600-5280 represents the wavelength range 
over which the slope is calculated. To avoid division by zero, implement

if meanS == 0:
    S = 0
    
"""

import os
import numpy as np
import pandas as pd
import astropy.units as u
from scipy.interpolate import interp1d
from mpdaf.obj import Cube, Spectrum

dates = ['2019-11-14', '2019-11-15', '2019-11-26', '2019-12-05', 
         '2019-12-06', '2019-12-21', '2019-12-23', '2019-12-29', 
         '2019-12-31', '2020-02-02', '2020-02-04', '2020-02-16', 
         '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']

def initialise(date):
    """
    Divides data cubes of 'date' by solar spectrum 'spectrum.dat'. Extracts
    three maps centered at 5280A,7080A,8600A of width 40A collapsed along 
    wavelength axis y taking the average. Co-adds maps with other exposures on
    same night by taking median. Saves the three maps belonging to date as
    .npy arrays and returns them as variables.

    Parameters
    ----------
    date : string

    Returns
    -------
    combinedmap1 : 2D numpy array
        blue map.
    combinedmap2 : 2D numpy array
        red map.
    combinedmap3 : 2D numpy array
        infrared map.

    """
    file = '../extra/spectrum.dat' #solar spectrum file
    df = pd.read_csv(file, names=['l', 'f', 'ignore'], delimiter=' ', skipinitialspace=(True), index_col=False) #l for lambda, f for flux
    
    bluemaps = []
    redmaps = []
    infraredmaps=[]
    directory = '../data/Shifted/' #driectory to shifted cubes
    for i, filename in enumerate(sorted(os.listdir(directory))):
        if filename[8:18] == date: #ensure date is in at indexed point in filename string
            cube = Cube(directory+filename)
            
            wavecoord = cube.wave.coord(unit=u.nm)
            fluxunit =  1e-20 * u.erg / (u.angstrom * u.s * u.cm **2)
            cube.unit = fluxunit
            
            df1 = df[(df.l >= 450) & (df.l <= 950)]
            solarflux = shift_wavelengths(df1.l, df1.f, wavecoord)
            fluxunit = u.W / ( u.nm * u.m ** 2)
            solarspec = Spectrum(wave=cube.wave, data=solarflux, unit=fluxunit)

            cube2 = cube / solarspec
            
            subcube1 = cube2.select_lambda(lbda_min=5260, lbda_max=5300, unit_wave=u.angstrom) #map 5260A-5300A
            map1 = subcube1.sum(axis=0) / 40
            image1 = map1.data
            bluemaps.append(image1)
            subcube2 = cube2.select_lambda(lbda_min=7060, lbda_max=7100, unit_wave=u.angstrom) #map 7060A-7100A
            map2 = subcube2.sum(axis=0) / 40
            image2 = map2.data
            redmaps.append(image2)
            subcube3 = cube2.select_lambda(lbda_min=8580, lbda_max=8620, unit_wave=u.angstrom) #map 8580A-8620A
            map3 = subcube3.sum(axis=0) / 40
            image3 = map3.data
            infraredmaps.append(image3)
            
    combinedmap1 = np.median(bluemaps, axis=0)
    combinedmap2 = np.median(redmaps, axis=0)
    combinedmap3 = np.median(infraredmaps, axis=0)
    
    np.save('../data/MapsForSlopes/'+'BlueMap_'+date+'.npy', combinedmap1)
    np.save('../data/MapsForSlopes/'+'RedMap_'+date+'.npy', combinedmap2)
    np.save('../data/MapsForSlopes/'+'InfraredMap_'+date+'.npy', combinedmap3)
    
    return combinedmap1, combinedmap2, combinedmap3

def shift_wavelengths(data_wavecoord, data, ideal_wavecoord):
    """
    Creates an interpolated spectrum from data_wavecoord and returns new 
    spectrum with wavecoord matching ideal_wavecoord.


    Parameters
    ----------
    data_wavecoord : 1D array
    data : 1D array
    ideal_wavecoord : 1D array
    Returns
    -------
    data_int_data : 1D array
        data of input variable data at new wavelength matching ideal_wavecoord

    """
    data_int = interp1d(data_wavecoord, data, kind='cubic',fill_value='extrapolate')        
    data_int_data = data_int(ideal_wavecoord)
    return data_int_data

def slope(bluemap, redmap, infraredmap, date):
    """
    Computes S' for each pixel and saves as .npy array
    
    S' = ((S2-S1)/(8600-5280)) * (1/meanS) * 100000
    

    Parameters
    ----------
    bluemap : 2D array
        S1 in above formula.
    redmap : 2D array
        meanS in above formula.
    infraredmap : 2D array
        S2 in above formula.
    date : string
        for saving file name.

    Returns
    -------
    None.

    """
    slopemap=np.empty((bluemap.shape[0], bluemap.shape[1]))
    for i in range(0,bluemap.shape[0]-1,1):
        #print('column: ',i,' being calculated')
        for j in range(0,bluemap.shape[1]-1,1):
            S1 = bluemap[i,j]
            S2 = infraredmap[i,j]
            meanS = redmap[i,j]
            
            if meanS == 0:
                S = 0
            else:
                S = ((S2-S1)/(8600-5280)) * (1/meanS) * 100000
            slopemap[i,j] = S
    np.save('../data/SlopeMaps/'+'ReflectivityGradients_'+date+'.npy', slopemap)
    return 

def main1():
    """
    Computes and saves numpy arrays of relectivity gradients from sctratch.

    Returns
    -------
    None.

    """
    for date in dates:
        bluemap, redmap, infraredmap = initialise(date)
        slope(bluemap, redmap, infraredmap, date)
    return

def main2():
    """
    Computes and saves numpy arrays of relectivity gradients from 
	already created and saved blue, red and infrared maps.

    Returns
    -------
    None.

    """
    for date in dates:
        bluemap = np.load('../data/MapsForSlopes/'+'BlueMap_'+date+'.npy')
        redmap = np.load('../data/MapsForSlopes/'+'RedMap_'+date+'.npy')
        infraredmap = np.load('../data/MapsForSlopes/'+'InfraredMap_'+date+'.npy')
        slope(bluemap, redmap, infraredmap, date)
    return

#main1()
#main2()


    
    
    
    
    
    
    
    
    
    
