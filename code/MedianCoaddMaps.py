#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 10:14:27 2023

@author: sed79

Combines maps by taking the median

Various date combinations can be used. Run function with info equal to one of 
the date combinations listed as global variable.

Move the created maps out of the GAS_ALL folders if they are a combination of dates.
"""

import os
import numpy as np
from astropy.io import fits

datesall= ['2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21', 
           '2019-12-23', '2019-12-29', '2019-12-31', '2020-02-02', '2020-02-04', 
           '2020-02-16', '2020-02-25', '2020-02-28', '2020-03-16', '2020-03-19']
presplit= ['2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21',
           '2019-12-23', '2019-12-29', '2019-12-31', '2020-02-02', '2020-02-04', 
           '2020-02-16', '2020-02-25', '2020-02-28']
postsplit = ['2020-03-16', '2020-03-19']
dates2019 = ['2019-11-15', '2019-11-26', '2019-12-05', '2019-12-06', '2019-12-21', 
           '2019-12-23', '2019-12-29', '2019-12-31']
dates2020 = ['2020-02-02', '2020-02-04', '2020-02-16', '2020-02-25', '2020-02-28', 
             '2020-03-16', '2020-03-19']
datesdec = ['2019-12-05', '2019-12-06', '2019-12-21', 
           '2019-12-23', '2019-12-29', '2019-12-31']
datesfeb = ['2020-02-02', '2020-02-04', '2020-02-16', '2020-02-25', '2020-02-28']
datesfinal7 = ['2020-02-25', '2020-02-28']

datesnov= ['2019-11-15', '2019-11-26']
datesearlydec = ['2019-12-05', '2019-12-06']


dates1 = (datesall, 'all.fits')
dates2 = (presplit, 'presplit.fits')
dates3 = (postsplit, 'postsplit.fits')
dates4 = (dates2019, '2019.fits')
dates5 = (dates2020, '2020.fits')
dates6 = (datesdec, 'dec.fits')
dates7 = (datesfeb, 'feb.fits')
dates8 = (datesfinal7, 'end_feb.fits')
dates9 = (datesnov, 'nov.fits')
dates99 = (datesearlydec, 'early_dec.fits')

dates10 = [['2019-11-15'], '2019-11-15.fits']
dates11 = [['2019-11-26'], '2019-11-26.fits']
dates12 = [['2019-12-05'], '2019-12-05.fits']
dates13 = [['2019-12-06'], '2019-12-06.fits']
dates14 = [['2019-12-21'], '2019-12-21.fits']
dates15 = [['2019-12-23'], '2019-12-23.fits']
dates16 = [['2019-12-29'], '2019-12-29.fits']
dates17 = [['2019-12-31'], '2019-12-31.fits']
dates18 = [['2020-02-02'], '2020-02-02.fits']
dates19 = [['2020-02-04'], '2020-02-04.fits']
dates20 = [['2020-02-16'], '2020-02-16.fits']
dates21 = [['2020-02-25'], '2020-02-25.fits']
dates22 = [['2020-02-28'], '2020-02-28.fits']
dates23 = [['2020-03-16'], '2020-03-16.fits']
dates24 = [['2020-03-19'], '2020-03-19.fits']
dates25 = [['2019-11-14'], '2019-11-14.fits']

def combine_CN(info):
    dates, datesname = info
    directory = '../data/ExtractedMaps/CN/'
    images = []
    for i, filename in enumerate(sorted(os.listdir(directory))):
        if filename[3:13] in dates:
            print('adding '+filename)
            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            image = hdulist[0].data
            images.append(image)
    combinedmap = np.median(images, axis=0)
    writing_name = 'CN_'+datesname
    hdu = fits.PrimaryHDU(combinedmap)
    hdu.writeto('../data/MedianCoaddedMaps/CN_ALL/'+writing_name)
    return
    
def combine_C2(info):
    dates, datesname = info
    directory = '../data/ExtractedMaps/C2/'
    images = []
    for i, filename in enumerate(sorted(os.listdir(directory))):
        if filename[3:13] in dates:
            print('adding '+filename)
            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            image = hdulist[0].data
            images.append(image)
    combinedmap = np.median(images, axis=0)
    writing_name = 'C2_'+datesname
    hdu = fits.PrimaryHDU(combinedmap)
    hdu.writeto('../data/MedianCoaddedMaps/C2_ALL/'+writing_name)
    return
    
def combine_NH2(info):
    dates, datesname = info
    directory = '../data/ExtractedMaps/NH2/'
    images = []
    for i, filename in enumerate(sorted(os.listdir(directory))):
        if filename[4:14] in dates:
            print('adding '+filename)
            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            image = hdulist[0].data
            images.append(image)
    combinedmap = np.median(images, axis=0)
    writing_name = 'NH2_'+datesname
    hdu = fits.PrimaryHDU(combinedmap)
    hdu.writeto('../data/MedianCoaddedMaps/NH2_ALL/'+writing_name)
    return

def combine_DUST(info):
    dates, datesname = info
    directory = '../data/ExtractedMaps/DUST/'
    images = []
    for i, filename in enumerate(sorted(os.listdir(directory))):
        if filename[5:15] in dates:
            print('adding '+filename)
            hdulist = fits.open(directory+filename,ignore_missing_end=True)
            image = hdulist[0].data
            images.append(image)
    combinedmap = np.median(images, axis=0)
    writing_name = 'DUST_'+datesname
    hdu = fits.PrimaryHDU(combinedmap)
    hdu.writeto('../data/MedianCoaddedMaps/DUST_ALL/'+writing_name)
    return