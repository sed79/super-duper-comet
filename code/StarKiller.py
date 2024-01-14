#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 20 10:55:43 2023

@author: sophie

reduce date by date with starkiller
"""

import os
from cuber import cuber

### Checking pysynphot is correctly confirgured
print(os.environ['PYSYN_CDBS'])
# if not enter into terminal: export PYSYN_CDBS=/home/sophie/Coding/CubeScene/PYSYN_CDBS_files/

path_cubes = '../data/CUBES/'

def reduce_cubes(year_month_day):
    """runs cuber() on files specfied by month and year"""
    for i, file in enumerate(os.listdir(path=path_cubes)):
        if file[18:28] == year_month_day: #choose month of cubes to work with
            print('Starting '+file)
            cuber(file=path_cubes+file,cal_maglim=20,run=True)
            
# date by date
#reduce_cubes('2020-03-19')
