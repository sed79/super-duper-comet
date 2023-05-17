#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 14:14:27 2023

@author: sed79

Make table for latex

Definition of info queried
https://astroquery.readthedocs.io/en/latest/api/astroquery.jplhorizons.HorizonsClass.html#astroquery.jplhorizons.HorizonsClass.ephemerides

To read extra queried items, such as elevation, go x =  JPLtableData['EL'] once run.

"""


import matplotlib.pyplot as plt
import numpy as np
import os
#---
from astropy.io import fits
from astropy.io import ascii
import astropy.time
from astropy.time import Time
import astropy.table 
from astropy.table import Table, Column, MaskedColumn
#---
from astroplan import Observer
from astroquery.jplhorizons import Horizons
#---
import pandas as pd
import csv

au = 149597870.7 #1au in km

# All dates

obsDate_lst = ['2019-11-14T08:03:10.611',
 '2019-11-14T08:19:45.880',
               '2019-11-14T08:31:36.260',
               '2019-11-14T08:48:23.334',
 '2019-11-15T08:02:02.188',
 '2019-11-15T08:18:33.659',
               '2019-11-15T08:30:28.508',
               '2019-11-15T08:53:20.124',
 '2019-11-26T07:17:24.150',
 '2019-11-26T07:34:33.386',
 '2019-11-26T07:46:28.548',
 '2019-11-26T08:03:33.586',
 '2019-12-05T07:09:26.388',
 '2019-12-05T07:26:38.436',
 '2019-12-05T07:38:31.418',
 '2019-12-05T07:55:01.516',
 '2019-12-06T07:37:32.770',
 '2019-12-06T07:53:57.483',
 '2019-12-06T08:05:52.444',
 '2019-12-06T08:25:06.054',
 '2019-12-21T07:56:19.257',
 '2019-12-21T08:13:29.848',
 '2019-12-21T08:25:25.648',
 '2019-12-21T08:42:26.639',
 '2019-12-23T05:50:43.552',
 '2019-12-23T06:07:20.214',
 '2019-12-23T06:19:16.900',
 '2019-12-23T06:42:59.231',
 '2019-12-29T05:18:41.077',
 '2019-12-29T05:38:15.405',
 '2019-12-29T05:50:13.126',
 '2019-12-29T06:09:06.900',
 '2019-12-31T05:36:51.322',
 '2019-12-31T05:53:46.784',
 '2019-12-31T06:05:41.099',
 '2019-12-31T06:22:04.886',
 '2020-02-02T04:17:53.363',
 '2020-02-02T04:34:27.258',
 '2020-02-02T04:46:27.238',
 '2020-02-02T05:03:48.161',
 '2020-02-04T03:18:11.907',
 '2020-02-04T03:34:48.409',
 '2020-02-04T03:46:44.140',
 '2020-02-04T04:05:56.771',
 '2020-02-16T04:07:13.749',
 '2020-02-16T04:24:06.115',
     '2020-02-16T04:35:59.862',
 '2020-02-16T04:51:50.792',
 '2020-02-16T05:04:46.042',
 '2020-02-25T03:28:18.901',
 '2020-02-25T03:43:47.143',
 '2020-02-25T03:55:37.851',
 '2020-02-25T04:11:53.150',
 '2020-02-28T04:02:19.394',
 '2020-02-28T04:20:27.879',
 '2020-02-28T04:32:18.164',
 '2020-02-28T04:48:45.942',
 '2020-03-16T04:41:49.600',
 '2020-03-16T04:58:16.132',
 '2020-03-16T05:10:08.903',
 '2020-03-16T05:28:34.842',
 '2020-03-19T05:23:42.010',
 '2020-03-19T05:46:49.016',
 '2020-03-19T05:58:44.909',
 '2020-03-19T06:14:45.202']



start_alt_angle = [28.743,
 32.394,
     34.978,
     38.591,
 29.29,
 32.936,
     35.542,
     40.467,
 28.025,
 31.922,
 34.623,
 38.484,
 33.244,
 37.156,
 39.861,
 43.617,
 40.4,
 44.137,
 46.851,
 51.221,
 54.763,
 58.552,
 61.183,
 64.934,
 28.822,
 32.347,
 34.896,
 39.976,
 26.183,
 30.151,
 32.589,
 36.463,
 31.102,
 34.507,
 36.916,
 40.237,
 33.612,
 35.928,
 37.591,
 39.965,
 26.173,
 28.45,
 30.091,
 32.726,
 36.912,
 38.805,
 40.104,
 41.786,
 43.109,
 35.373,
 36.937,
 38.106,
 39.662,
 39.441,
 41.03,
 42.015,
 43.309,
 44.268,
 44.832,
 45.159,
 45.531,
 45.278,
 45.398,
 45.353,
 45.183] #sourced from FITS headers

DIMM = [0.8,
0.8,
0.72,
0.72,
0.87,
0.68,
0.76,
0.57,
0.43,
0.74,
0.8,
0.74,
0.93,
1.34,
1.17,
1.21,
1.04,
1.61,
1.35,
1.47,
0.49,
0.57,
0.65,
0.51,
0.39,
0.28,
0.44,
0.32,
0.66,
0.71,
0.6,
0.48,
0.41,
0.38,
0.41,
0.32,
0.45,
0.53,
0.64,
0.54,
0.55,
0.66,
0.44,
0.45,
0.6,
0.79,
0.66,
0.47,
0.47,
0.81,
0.61,
0.75,
0.49,
0.72,
0.65,
0.6,
1.22,
1.23,
1.04,
1.09,
1.11,
0.61,
0.32,
0.36,
0.43] # sourced from night reports

# convert epoch date from 'fits' type to julian date
jd_obsDate_lst = []
for i in obsDate_lst:
    t = Time(i, format='fits', scale='utc')
    jd_obsDate_lst.append(t.jd) 
# query information from JPL Horizons database
observations = Horizons(id='C/2019 Q4', location='809',epochs=jd_obsDate_lst)
jpl_data = observations.ephemerides()
    
# JPL data 
JPLtableData = jpl_data['datetime_str','datetime_jd','RA','DEC','RA_rate','DEC_rate','Tmag','EL', 'airmass','r','r_rate','delta','delta_rate','elong','alpha', 'velocityPA', 'sunTargetPA']
JPL_colnames = ['datetime_str','datetime_jd','RA','DEC','RA_rate','DEC_rate','Tmag','elevation','airmass','r','r_rate','delta','delta_rate','elong','alpha', 'velocityPA', 'sunTargetPA']
JPL_astropy_tbl = Table(JPLtableData,names=JPL_colnames)

#compute km of 1" for each obs date
scale = JPLtableData['delta'] * np.tan(np.pi / (180*3600)) * au


# FITS data
FITStabledata = [start_alt_angle, DIMM, scale]
FITS_colnames = ['start_alt_angle', 'DIMM', 'scale']
FITS_astropy_tbl = Table(FITStabledata,names=FITS_colnames)
# combine jpl data and fits data
JPL_FITS_tbl = astropy.table.hstack((JPL_astropy_tbl,FITS_astropy_tbl))

# change name
final_combined_table = JPL_FITS_tbl

# we want the datetimes to appear nicely in the Latex table output
datetime_str_nice = []
for i in final_combined_table['datetime_str']:
    datetime_str_nice.append(i[0:-4]) #removing the last 4 characters of the datetime


#making a Night Reports table which includes the nice datetime
latex_tbl = Table([datetime_str_nice],names=["datetime_str_nice"])
latex_new_tbl = astropy.table.hstack((latex_tbl,JPL_FITS_tbl))
#joining this with the final table
latex_final_combined_table = latex_new_tbl

#WRITING TABLES

#as a CSV
ascii.write(final_combined_table,"simple_table_test.csv",format="csv",overwrite=True)
# just exporting scale column
ascii.write(latex_final_combined_table['datetime_str_nice','start_alt_angle','airmass','DIMM','r','delta','elong','alpha', 'sunTargetPA', 'velocityPA', 'scale'],"ObservationTable.csv",format="csv",overwrite=True,
            formats={'start_alt_angle':'%12.1f',
                     'airmass':'%12.2f',
                     'DIMM':'%12.2f',
                     'r':'%12.2f',
                     'delta':'%12.2f',
                     'elong':'%12.2f',
                     'alpha':'%12.2f',
                    'sunTargetPA':'%12.2f',
                    'velocityPA':'%12.2f',
                    'scale':'%12.0f'})

#as Latex code
ascii.write(latex_final_combined_table['datetime_str_nice','start_alt_angle','airmass','DIMM','r','delta','elong','alpha', 'sunTargetPA', 'velocityPA', 'scale'],format="latex",
           formats={'start_alt_angle':'%12.1f',
                    'airmass':'%12.2f',
                    'DIMM':'%12.2f',
                    'r':'%12.2f',
                    'delta':'%12.2f',
                    'elong':'%12.2f',
                    'alpha':'%12.2f',
                   'sunTargetPA':'%12.2f',
                   'velocityPA':'%12.2f',
                   'scale':'%12.0f'})