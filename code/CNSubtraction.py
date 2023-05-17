#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 11:00:54 2023

@author: sed79

This will take a dust subtracted data cube and further subtract the deep 
feature over the CN region, producing a subcube covering 8800A-9300A

Method to remove feature is;
Fitting a polynomial to 5 pixel radius aperture spectrum of 2019-12-05 07:09:26.388.
Using this polynomial as a model continuum
Scaling this continuum to each spaxel to then subtract
Scaling is done by fitting continuum to small range of wavelengths in spectrum that avoids 
the actual CN emission
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from mpdaf.obj import Cube, Spectrum
import astropy.units as u
from scipy.interpolate import interp1d

directory = '../data/DustSubtracted/'# Directory to dust subtracted cubes

#lists to divide up the processing
list1 = [0, 6, 12]
list2 = [1, 7, 13]
list3 = [2, 8, 14]
list4 = [3, 9, 15]
list5 = [4, 10, 16]
list6 = [5, 11, 17]


def main():
    modelcont = gen_model_cont()
    for i, filename in enumerate(sorted(os.listdir(directory))):
        if i in list1:
            print('beginning subtraction of '+filename)
            sub_save_cube(modelcont, filename)
    return

def gen_model_cont():
    """
    Creates model continuum object covering 8800A-9300A

    Returns
    -------
    cont : obj.spectrum.Spectrum

    """
    #get aperture spectrum and trim to 8800A-9300A
    cube = Cube('../data/DustSubtracted/DustSubtracted_2019-12-05T07:09:26.388.fits')
    specfull = cube.aperture((200.5,200.5), 5, unit_center=None, unit_radius=None, is_sum=False)
    specnogap = specfull.subspec(lmin=8800.00, lmax=9300) #specfull.get_end())
    
    #make gap spec (called spec) covering 8800A-9100A and 9230A-9300A
    spec = specnogap.copy()
    spec.mask_region(lmin=9100, lmax=9230, inside=True)
    
    #fit polynomial and create continuum spectrum
    z = spec.poly_fit(7)
    cont = spec.clone()
    cont.poly_val(z) 
    
    #plot the spectrum, continuum, and CN region
    fig, ax = plt.subplots()
    spec.plot(ax = ax)
    cont.plot(ax = ax)
    specnogap.plot(ax = ax, lmin=9135, lmax=9210, c='k')
    plt.grid()
    plt.title('2019-12-05')
    return cont

def shift_wavelengths(spectoshift, spectocopy):
    """
    Creates an interpolated spectrum from spectoshift and returns new spectrum with wavecoord matching spectocopy

    Parameters
    ----------
    spectoshift : obj.spectrum.Spectrum
    spectocopy : obj.spectrum.Spectrum

    Returns
    -------
    int_spectrum : obj.spectrum.Spectrum

    """
    spec_int = interp1d(spectoshift.wave.coord(), spectoshift.data, kind='cubic',fill_value='extrapolate')        
    spec_int_flux = spec_int(spectocopy.wave.coord())
    int_spectrum = Spectrum(wave=spectocopy.wave, data=spec_int_flux)
    return int_spectrum

def sub_save_cube(cont, filename):
    """
    Creates new data cube with model continuum subtracted from CN region. 
    Extracts map from cube.

    Parameters
    ----------
    cont : obj.spectrum.Spectrum
        The model continuum to subtract.

    Returns
    -------
    None.

    """

    #establish cube and empty subcube
    cube = Cube(directory+filename)
    subcube = cube.select_lambda(lbda_min=8800, lbda_max=9300, unit_wave=u.angstrom)
    reconstructed_data = np.empty((subcube.shape[0],subcube.shape[1],subcube.shape[2])) #empty cube over wavelengths 8800A-9300A
    continuum = shift_wavelengths(cont, subcube[:,0,0]) #interpolate model continuum to match cube's wavecoord
    
    #iterate through cube, subtract, and fill empty subcube
    for i in range(0,subcube.shape[1]-1,1):
        print('column: ',i,' being subtracted')
        for j in range(0,subcube.shape[2]-1,1):
            #get sub spectra
            spec = subcube[:, i, j] #specific spaxel
            specside = spec.subspec(lmin=9080.00, lmax=9100)
            specCN = spec.subspec(lmin=9135.00, lmax=9200.00)
            mean = np.mean(specside.data)
            #check if cont needs resampling
            if np.array_equal(spec.wave.coord(), continuum.wave.coord()) == False:
                continuum = shift_wavelengths(cont, spec) #interp cont to match cube's wavecoord
            #scale and sub
            y = continuum.subspec(lmin=9090) #midpoint of sidewall
            scale = mean / y
            contscaled = scale * continuum
            subbed = spec - contscaled
            #fill empty subcube
            reconstructed_data[:,i,j] = subbed.data
            
    
     
    #save subtracted subcube
    newcube = Cube(data=reconstructed_data,primary_header=cube.primary_header,data_header=cube.data_header,wave=subcube.wave,wcs=subcube.wcs)
    newcube.write('../data/CNSubtracted/CNSub_'+filename[15:])
    return

main() 
