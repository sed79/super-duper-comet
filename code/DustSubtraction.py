#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 10:08:05 2023

@author: sed79

Completes the dust subtraction of each cube. Cubes have previously been clipped of extreme pixels and had their centre shifted.
"""

import os
import numpy as np
from scipy.interpolate import interp1d
from astropy.modeling import models, fitting
from mpdaf.obj import Cube

directory = '../data/Shifted/' #direct for shifted cubes
dust_txt = '../extra/67PSpectrum_Dop.txt' #dpectrum of 67P

#lists to divide up the processing
list1 = [3, 29]
list2 = [4, 38]
list3 = [8, 39]
list4 = [9]
list5 = [13]
list6 = [14]
list7 = [23]
list8 = [24]
list9 = [28]

def main():
    """
    Makes new, dust-subtracted, .fits files for each of the .fits files in 
    the directory. 

    """
    dust = np.loadtxt(dust_txt)
    dust_spectrum = (dust[:,0], dust[:,1])
    #iterate through directory of cubes
    for i, filename in enumerate(sorted(os.listdir(directory))):
        file = directory + filename
        if i in list1: #dividing up work load
            scidata, comet_wavecoord = get_data(file)
            dust_interpolated_spectrum, dust_int = shift_wavelengths(dust_spectrum, comet_wavecoord)
            dust_gas_free_spectrum = remove_gas_emission(dust_interpolated_spectrum)
            reconstructed_data = subtract_dust(scidata, comet_wavecoord, dust_interpolated_spectrum, dust_gas_free_spectrum)
            new_fits(reconstructed_data, filename)
    return
            
def get_data(file):  
    """
    Extracts science data (scidata) and spectral coordinate axis (wavecoord) 
    from the specified .fits file

    Parameters
    ----------
    filename : string. Name of file in directory to extract data from. 

    Returns
    -------
    scidata : 3D numpy array. The 3D cube of data from the .fits file. 
    comet_wave_axis : 1D numpy array. The wavelengths corresponding to the 
                      observed data in the .fits file.

    """
    cube = Cube(file)
    scidata = cube.data
    wavecoord = cube.wave.coord()
    return scidata, wavecoord

def shift_wavelengths(dust_spectrum, comet_wavelengths):
    """
    Creates an interpolated dust spectrum so that a new dust spectrum, with 
    values measured at the same wavelengths as the observed comet, can be 
    created.

    Parameters
    ----------
    dust_spectrum : tuple. First entry is 1D array of wavelengths the comet 67P 
                    was observed at and the second entry is a 1D array the 
                    corresponding flux of 67P.
    comet_wavelengths : 1D numpy array. Wavelengths coresponding to observations 
                        of 2I for the specific observation being dust subtracted.

    Returns
    -------
    dust_int_spectrum : tuple. First entry is 1D array of wavelengths 2I 
                        was observed for the specific observation being dust 
                        subtracted at and the second entry is a 1D array the 
                        corresponding flux of 67P.

    """
    dust_int_wavelengths = np.copy(comet_wavelengths)
    dust_int = interp1d(dust_spectrum[0], dust_spectrum[1], kind='cubic',fill_value='extrapolate')        
    dust_int_flux = dust_int(comet_wavelengths)
    dust_int_spectrum = (dust_int_wavelengths, dust_int_flux)
    return dust_int_spectrum, dust_int

def remove_gas_emission(spectrum):
    """
    Removes the wavelengths and flux values of the input spectrum corresponding 
    to gas emisison to return a 'gas-free' spectrum. The gas-free wavelength
    ranges (in angstroms) are 4780.0-4830.0, 5250.0-5300.0, 
    5800.0-5850.0, 6400.0-6450.0, 6765.0-6815.0, 7130.0-7180.0, 
    and 8560.0-8610.0.

    Parameters
    ----------
    spectrum : tuple. First entry is 1D array of wavelengths and the second 
                      entry is a 1D array the flux.

    Returns
    -------
    gas_free_spectrum : tuple. First entry is 1D array of gass-free wavelengths 
                        and the second entry is a 1D array the corresponding flux.

    """
    wavelengths = spectrum[0]
    flux = spectrum[1]
    gas_free_wavelengths = np.asarray([])
    gas_free_flux = np.asarray([])
    
    for i in range(len(wavelengths)):
        if 4780.0<wavelengths[i]<4830.0 or 5250.0<wavelengths[i]<5300.0 or 5800.0<wavelengths[i]<5850.0  or 6400.0<wavelengths[i]<6450.0 or 6765.0<wavelengths[i]<6815.0 or  7130.0<wavelengths[i]<7180.0 or 8560.0<wavelengths[i]<8610.0 or 9010.0<wavelengths[i]<9030.0:
            gas_free_wavelengths = np.append(gas_free_wavelengths,wavelengths[i])
            gas_free_flux = np.append(gas_free_flux,flux[i])
    gas_free_spectrum = (gas_free_wavelengths, gas_free_flux)
    return gas_free_spectrum

def make_model(comet_spectrum, dust_spectrum):
    """
    Creates a function which models the shape of the flux of the gas-free 
    comet emission divided by the gas-free dust emisison. The function is a
    2nd order polynomial that was fitted using the Levenberg-Marquardt 
    algorithm and least squares statistics.

    Parameters
    ----------
    comet_spectrum : tuple. First entry is 1D array of gas-free wavelengths 
                     and the second entry is a 1D array of the corresponding 
                     comet flux.
    dust_spectrum : tuple. First entry is 1D array of gas-free wavelengths 
                    and the second entry is a 1D array the corresponding 
                    dust flux.

    Returns
    -------
    gDust : 1D polynomial model. 

    """
    wavelengths = comet_spectrum[0]
    comet_flux = comet_spectrum[1]
    dust_flux = dust_spectrum[1]
    div_flux = comet_flux/dust_flux
    
    g_initDust = models.Polynomial1D(2)
    fit_gDust = fitting.LinearLSQFitter()
    gDust = fit_gDust(g_initDust, wavelengths, div_flux)
    return gDust, div_flux
    
def correct_slope(spectrum, model):
    """
    Corrects the slope of the spectrum by multiplying the flux by 
    the model. Model takes input of the same wavelengths as the spectrum 
    to produce a 1D array of the polynomial before multiplication.

    Parameters
    ----------
    spectrum : tuple. First entry is 1D array of wavelengths and the second 
               entry is a 1D array of the corresponding flux.
    model : 1D polynomial model.

    Returns
    -------
    slope_corrected_spectrum : tuple. First entry is 1D array of wavelengths 
                                and the second entry is a 1D array of the 
                                corresponding slope-corrected flux.

    """
    wavelengths = spectrum[0]
    flux = spectrum[1]
    fitted_line = np.zeros_like(wavelengths)
    slope_corrected_flux= np.zeros_like(flux)

    fitted_line = model(wavelengths)
    slope_corrected_flux = flux * fitted_line
    
    slope_corrected_spectrum = (wavelengths, slope_corrected_flux)
    return slope_corrected_spectrum
    
def subtract_dust(scidata, comet_wavelengths, dust_interpolated_spectrum, dust_gas_free_spectrum):
    """
    Subtracts the dust from the flux in the cube. This is done by:
    
    Iterating over each spaxel, normalising the dust flux to the spaxel's flux 
    by multiplcation of relative total gas-free emission of 2I to 67P, 
    subtracting the normalised dust flux from the comet's flux, and making
    a new 3D array of the subtracted spectra.

    Parameters
    ----------
    scidata : 3D numpy array. The 3D cube of data from the .fits file. 
    slope_corrected_dust_spectrum : tuple. First entry is 1D array of wavelengths 
                                and the second entry is a 1D array of the 
                                corresponding slope-corrected dust flux.

    Returns
    -------
    reconstructed_data : 3D numpy array. The cube of data with dust-subtacted flux.
    
    """
    wavelengths = comet_wavelengths
    reconstructed_data = np.empty((scidata.shape[0],scidata.shape[1],scidata.shape[2]))
    
    for i in range(0,scidata.shape[1]-1,1):
        print('column: ',i,' being subtracted')
        for j in range(0,scidata.shape[2]-1,1):
                
                flux = scidata[:,i,j]
        
                comet_gas_free_spectrum = remove_gas_emission((wavelengths, flux))
                    
                model, div_flux = make_model(comet_gas_free_spectrum, dust_gas_free_spectrum)
                slope_corrected_dust_spectrum = correct_slope(dust_interpolated_spectrum, model)
                gas_free_slope_corrected_dust_spectrum = remove_gas_emission(slope_corrected_dust_spectrum)
            
                dust_flux = slope_corrected_dust_spectrum[1]
                summed_dust_flux = np.sum(gas_free_slope_corrected_dust_spectrum[1])
                summed_comet_flux = np.sum(comet_gas_free_spectrum[1])
                if summed_comet_flux != 0:
                    normalised_dust = (dust_flux / summed_dust_flux) * summed_comet_flux
                else:
                    normalised_dust = dust_flux * 0
                    
                subtracted_data = flux - normalised_dust
                reconstructed_data[:,i,j] = subtracted_data
                    
    return reconstructed_data


def new_fits(reconstructed_data, filename):
    """
    Writes the reconstructed data to a new .fits file into the working 
    directory with the same header information as the file it was made from and 
    new filename 'DustSubtracted_existingfilename'.

    Parameters
    ----------
    reconstructed_data : 3D numpy array.
    filename : string. 

    Returns
    -------
    None.

    """
    cube=Cube(directory+filename)
    cube_data_only = Cube(data=reconstructed_data,primary_header=cube.primary_header,data_header=cube.data_header,wave=cube.wave,wcs=cube.wcs)
    cube_data_only.write('../data/DustSubtracted/DustSubtracted_'+filename[8:]) #ensure filename[8:] is a good name to save with

main()
