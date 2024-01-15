#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 14:23:47 2023

@author: sed79

Will shift cubes to put comet centre at 200,200. Cubes must be in list 
2I_Center_List.txt or 2I_Center_List_Extra.txt

"""

import numpy as np
import astropy.units as u
from mpdaf.obj import Cube
from mpdaf.obj import iter_ima

directory = '../data/ExtremesMasked/' #directory of cubes to be shifted
center_list = '../extra/2I_Center_List_master_updated_diffcubes.txt' # cubes to shift

#lists to divide up the processing
list1 = [0, 1, 2, 3, 4]
list2 = [5, 6, 7, 8, 9]
list3 = [10, 11, 12, 13, 14]
list4 = [15, 16, 17, 18, 19]
list5 = [20, 21, 22, 23, 24]
list6 = [25, 26, 27, 28, 29]
list7 = [30, 31, 32, 33, 34]
list8 = [35, 36, 37, 38, 39]
list9 = [40, 41, 42]

def shift_all():
    """
    Shifts each cube listed in 2I_Center_List.txt using shift_start.

    Returns
    -------
    None.

    """
    names = np.genfromtxt(center_list,usecols=0,dtype=str,encoding='utf-8-sig')
    xcen = np.genfromtxt(center_list,usecols=1,dtype=float)
    ycen = np.genfromtxt(center_list,usecols=2,dtype=float)    

    for i, name in enumerate(names):
        if i in list1: #diving up work load
            print('starting: '+name)
            comet_center = (xcen[i], ycen[i])
            shift_start(name, comet_center)
            print('finished: '+name)
    return
        
def shift_start(name, comet_center,k_dbg=None):
    """
    Shifts the comet center onto the center of the pixel 200, 200. 
    
    Sets values of outer spatial borders with thickness 15 pix to 0.0. This 
    removes the large borders of masked values. The shifting then fills the 
    remaining masked data with a spline interpolation between neighbouring 
    values. Order of this interpolation is not able to be specified, and 
    documentation does not say. JK linear is faster/works.
    
    Creates empty data cube and places shifted cube into it so that comet center 
    is in middle of 200,200 spaxel. Cube is shifted using 3rd order spline 
    interpolation.Fills rest of new data cube spatial dimensions 
    external to original spatial 
    dimensions (eg the area outside of the (307, 319) within the (400,400)) 
    with empty array values.
    
    Writes new .fits file named 'Shifted_dateobserved.fits', which contains 
    cube of shifted data with dimensions 400,400. 
    
    Ensure cube's rotation of y axis wrt DEC is -0.0. 

    Parameters
    ----------
    name : string
    comet_center : tuple

    Returns
    -------
    None.

    """
    xcen, ycen = comet_center
    cube = Cube(directory+'Masked_'+name[18:])
    writename = 'Shifted_'+name[18:]
    
    start_y = cube.get_start(unit_wave=u.nm, unit_wcs=u.deg)[1] #Return [lbda,y,x] at the center of pixel (0,0,0)
    start_x = cube.get_start(unit_wave=u.nm, unit_wcs=u.deg)[2]
    step_y = cube.get_step(unit_wave=u.nm, unit_wcs=u.deg)[1] #Return the cube steps [dlbda,dy,dx]
    width = cube.get_range(unit_wave=u.nm, unit_wcs=u.deg)[-1] - cube.get_range(unit_wave=u.nm, unit_wcs=u.deg)[2]
    xlen = cube.shape[2]
    step_x = -1*width/(xlen-1)
    ycen_deg = start_y+(step_y*(ycen-0.5))
    xcen_deg = start_x+(step_x*(xcen-0.5))
    
    newstart_y = ycen_deg - (200*step_y)
    newstart_x = xcen_deg - (200*step_x)
    newstart = (newstart_y, newstart_x)
    newdim=(400, 400) #for unchanged shape, newdim = (cube.shape[1], cube.shape[2])
    newstep = (cube.get_step(unit_wave=u.nm, unit_wcs=u.deg)[1], cube.get_step(unit_wave=u.nm, unit_wcs=u.deg)[2])
    
    # set bordering values to 0.0, unmasking the values
    cube.data[:, :15, :]=0.0
    cube.data[:, -15:, :]=0.0
    cube.data[:, :, :15]=0.0
    cube.data[:, :, -15:]=0.0
    
    cube2 = Cube.new_from_obj(cube, data=np.empty((cube.shape[0], 400, 400)), var=False)
    print('Starting resampling')

    for ima, k in iter_ima(cube, index=True): 
        if k%100 == 0:
            print('Shifted '+str(k)+' images')
        cube2[k,:,:] = ima.resample(newdim, newstart, newstep, order=3, unit_step=u.deg, interp='linear')
    cube2.write('../data/Shifted/'+writename)
    return
    
shift_all()
