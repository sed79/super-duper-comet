#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 11:24:02 2023

@author: sed79

Makes a trajoctory plot of 2I with data points corresponding to dates when 2I was observed by MUSE.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.path as mpath
np.set_printoptions(precision=3)
from jplephem.spk import SPK
from spktype21 import SPKType21

AU = 149597871 #1 AU in km

# all observations
table = pd.read_csv('../extra/JPL_query_table.csv', delimiter = ',')
jds_all = table['datetime_jd']

# data subset - one filename per night 
subset_names = np.genfromtxt('../extra/2I_one_jd_of_each_obs_epoch.txt',usecols=0,dtype=str,encoding='utf-8-sig')
for i in range(len(subset_names)):
    subset_names[i] = subset_names[i][18:41]
subset_table = table[table['obs_date'].isin(subset_names)]
jds_subset = subset_table['datetime_jd']

# Special dates
jds_disc = [2458725.50000] #JD of discovery at beginning of Aug 30
jds_peri =  [2458823.850764514] #jd of perihelion
jds_ob = [2458915.50000] #jd of outburst mar 6-9 2019, chosen mar 7
jds_split = [2458938.50000] #jd resolving of splitting event march 30 2019
jds_year = np.linspace(2458484.50000, 2459215.50000, 100) #the years 2019-2021
jds_month = np.linspace(2458484.50000, 2458604.500002, 100) #the months 01-2019-05-2019
jds_2I = np.linspace(2458710.50000, 2458939.50, 365) #2019-08-15 to  (querey limits are Start=2019-08-01 TDB , Stop=2020-04-01)

#%% Generate Positions

kernel1 = SPKType21.open('../extra/1003639 (3).bsp') #2I
kernel2 = SPK.open('../extra/de430.bsp') #solar system

x1 = [] #2I wrt Solar System Barycenter (icecream data)
y1 = []
z1 = []
x2 = [] #2I wrt Solar System Barycenter (gravel data)
y2 = []
z2 = []
x3 = [] #2I wrt Solar System Barycenter (almost year long data)
y3 = []
z3 = []
x4 = [] #Sun wrt Solar System Barycenter (year long data)
y4 = []
z4 = []
x5 = [] #Earth wrt Solar System Barycenter (year long data)
y5 = []
z5 = []
x6 = [] #Mars wrt Solar System Barycenter (year long data)
y6 = []
z6 = []
x7 = [] #Jupiter wrt Solar System Barycenter (year long data)
y7 = []
z7 = []
x8 = [] #Mercury wrt Solar System Barycenter (year long data)
y8 = []
z8 = []
x9 = [] #Venus wrt Solar System Barycenter (year long data)
y9 = []
z9 = []
x10 = [] #2I wrt Solar System Barycenter (at discovery)
y10 = []
z10 = []
x110 = [] #Earth wrt Solar System Barycenter (over 2I's observations ice cream)
y110 = []
z110 = []
x111 = [] #Earth wrt Solar System Barycenter (over 2I's observations gravel)
y111 = []
z111 = []
x13 = [] #Earth wrt Solar System Barycenter (at discovery)
y13 = []
z13 = []
x14 = [] #2I wrt Solar System Barycenter (at splitting event)
y14 = []
z14 = []
x15 = [] #Earth wrt Solar System Barycenter (at splitting event)
y15 = []
z15 = []
x16 = [] #2I wrt Solar System Barycenter (at outburst event)
y16 = []
z16 = []
x17 = [] #Earth wrt Solar System Barycenter (at outburst event)
y17 = []
z17 = []

for jd in jds_subset:
    positionA, velocityA = kernel1.compute_type21(10,1003639, jd) # 2I wrt sun
    positionB = kernel2[0, 10].compute(jd) # sun wrt ssb 
    position1 = (positionA - positionB) #2I wrt ssb 
    position110 = kernel2[0, 3].compute(jd) #Earth wrt ssb
    x1.append(position1[0]/AU)
    y1.append(position1[1]/AU)
    z1.append(position1[2]/AU)
    x110.append(position110[0]/AU)
    y110.append(position110[1]/AU)
    z110.append(position110[2]/AU)
    

for jd in jds_2I:
    positionA, velocityA = kernel1.compute_type21(10,1003639, jd) # 2I wrt sun
    positionB = kernel2[0, 10].compute(jd) # sun wrt ssb 
    position3 = (positionA - positionB) #2I wrt ssb
    x3.append(position3[0]/AU)
    y3.append(position3[1]/AU)
    z3.append(position3[2]/AU)

for jd in jds_year:
    position4 = kernel2[0, 10].compute(jd) # sun wrt ssb 
    position5 = kernel2[0, 3].compute(jd) # earth wrt ssb
    position6 = kernel2[0, 4].compute(jd) # mars wrt ssb
    position7 = kernel2[0, 5].compute(jd) # jupiter wrt ssb
    
    position9 = kernel2[0, 2].compute(jd) # venus wrt ssb
    x4.append(position4[0]/AU)
    y4.append(position4[1]/AU)
    z4.append(position4[2]/AU)
    x5.append(position5[0]/AU)
    y5.append(position5[1]/AU)
    z5.append(position5[2]/AU)
    x6.append(position6[0]/AU)
    y6.append(position6[1]/AU)
    z6.append(position6[2]/AU)
    x7.append(position7[0]/AU)
    y7.append(position7[1]/AU)
    z7.append(position7[2]/AU)
    
    x9.append(position9[0]/AU)
    y9.append(position9[1]/AU)
    z9.append(position9[2]/AU)
    
for jd in jds_month:
    position8 = kernel2[0, 1].compute(jd) # mercury wrt ssb
    x8.append(position8[0]/AU)
    y8.append(position8[1]/AU)
    z8.append(position8[2]/AU)
    
for jd in jds_disc:
    positionA, velocityA = kernel1.compute_type21(10,1003639, jd) # 2I wrt sun
    positionB = kernel2[0, 10].compute(jd) # sun wrt ssb 
    position10 = (positionA - positionB) #2I wrt ssb
    position13 = kernel2[0, 3].compute(jd) # earth wrt ssb
    x10.append(position10[0]/AU)
    y10.append(position10[1]/AU)
    z10.append(position10[2]/AU)
    x13.append(position13[0]/AU)
    y13.append(position13[1]/AU)
    z13.append(position13[2]/AU)
    
for jd in jds_split:
    positionA, velocityA = kernel1.compute_type21(10,1003639, jd) # 2I wrt sun
    positionB = kernel2[0, 10].compute(jd) # sun wrt ssb 
    position14 = (positionA - positionB) #2I wrt ssb
    position15 = kernel2[0, 3].compute(jd) # earth wrt ssb
    x14.append(position14[0]/AU)
    y14.append(position14[1]/AU)
    z14.append(position14[2]/AU)
    x15.append(position15[0]/AU)
    y15.append(position15[1]/AU)
    z15.append(position15[2]/AU)
  
for jd in jds_ob:
    positionA, velocityA = kernel1.compute_type21(10,1003639, jd) # 2I wrt sun
    positionB = kernel2[0, 10].compute(jd) # sun wrt ssb 
    position16 = (positionA - positionB) #2I wrt ssb
    position17 = kernel2[0, 3].compute(jd) # earth wrt ssb
    x16.append(position16[0]/AU)
    y16.append(position16[1]/AU)
    z16.append(position16[2]/AU)
    x17.append(position17[0]/AU)
    y17.append(position17[1]/AU)
    z17.append(position17[2]/AU)

kernel1.close()
kernel2.close()

#%% Plot 

def plot(gravel=False, save=False, withoutbackground=False):
    from mpl_toolkits.mplot3d.axes3d import Axes3D
    
    fig = plt.figure(figsize=(16,16))
    ax = fig.add_subplot(projection='3d', computed_zorder=False)
    ax.set_proj_type('persp')  # FOV = 90 deg
    
    # Subset data
    for i in range(len(x1)):
        ax.scatter3D(x1[i], y1[i], z1[i], c='olivedrab', edgecolors=('k'), marker='o', s=60.0, linewidth=1, label='Obs Epochs', zorder=6, depthshade = False)
    ax.scatter3D(x110, y110, z110, c='olivedrab', marker='o', s=60.0, edgecolor='k', linewidth=1, label='Earth at obs epochs', zorder=6, depthshade = False)
        
    # 2I at special dates
    ax.scatter3D(x10, y10, z10, c='w', marker='D', s=60.0, edgecolor='k', linewidth=1, label='Discovery', zorder=5, depthshade = False)
    ax.scatter3D(x13, y13, z13, c='w', marker='D', s=60.0, edgecolor='k', linewidth=1, label='Discovery', zorder=5, depthshade = False)
    ax.scatter3D(x14, y14, z14, c='w', marker='s', s=60.0, edgecolor='k', linewidth=1, label='Splitting', zorder=5, depthshade = False)
    ax.scatter3D(x15, y15, z15, c='w', marker='s', s=60.0, edgecolor='k', linewidth=1, label='Splitting', zorder=5, depthshade = False)
    ax.scatter3D(x16, y16, z16, c='w', marker='^', s=70.0, edgecolor='k', linewidth=1, label='Outburst', zorder=7, depthshade = False)
    ax.scatter3D(x17, y17, z17, c='w', marker='^', s=70.0, edgecolor='k', linewidth=1, label='Outburst', zorder=7, depthshade = False)
    
    # 2I's trajectory
    ax.plot3D(x3, y3, z3, label = '2I path', c='grey', linestyle='dashed', linewidth=2)
    
    # Planet Rings
    ax.scatter3D(x4, y4, z4, c='gold',  label='Sun', edgecolors=['xkcd:brick orange'], linewidth=0.5, s=200.0)
    ax.plot3D(x8, y8, z8, c='k', label='Mercury')
    ax.plot3D(x9, y9, z9, c='k', label='Venus')
    ax.plot3D(x5, y5, z5, c='k', label='Earth')
    ax.plot3D(x6, y6, z6, c='k', label='Mars')
    #ax.plot3D(x7, y7, z7, c='orange', label='Jupiter')
    
    # First remove fill
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    
    # Now set color to white (or whatever is "invisible")
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')
    
    # Bonus: To get rid of the grid as well:
    ax.xaxis._axinfo["grid"].update({"linewidth":0.1, 'color':'gainsboro'})
    ax.yaxis._axinfo["grid"].update({"linewidth":0.1,  'color':'gainsboro'})
    ax.zaxis._axinfo["grid"].update({"linewidth":0.1, 'color':'gainsboro'})
    
    
    # Axis labels
    ax.set_xlabel('(au)', labelpad = 5, fontsize='12', rotation=0)
    ax.set_ylabel('(au)', labelpad = 7, fontsize='12', rotation=0)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel("(au)", labelpad = 10, fontsize='12', rotation=90)
    ax.set(xlim=(-2.0, 1.5), ylim=(-2.0, 2.0))
    ax.set(yticklabels=(-2.0, '', -1.0, '', 0.0, '', 1.0, '', 2.0), xticklabels=(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0), zticklabels=('', -2.0, -1.0, 0.0, 1.0, 2.0))
    
    # Scale image
    x_scale=1
    y_scale=1
    z_scale=1 #1.5
    scale=np.diag([x_scale, y_scale, z_scale, 1.0])
    scale=scale*(1.0/scale.max())
    scale[3,3]=1.0
    def short_proj():
      return np.dot(Axes3D.get_proj(ax), scale)
    ax.get_proj=short_proj
    
    
    #View angle
    ax.view_init(5,235) #tilt up-down, spin around, best(5, 235)
    
    from mpl_toolkits.mplot3d import axes3d
    import mpl_toolkits.mplot3d.art3d as art3d
    from matplotlib.patches import Arrow, PathPatch
    arrowfig = Arrow(-2,2.2,-0.08,-0.5,width =0.1, color = 'k', alpha=0.8) #x,z,dx,dz
    ax.add_patch(arrowfig)
    art3d.pathpatch_2d_to_3d(arrowfig, z=0.6, zdir="y")
    
    arrowfig = Arrow(-1.8,-2.2,0.15,-0.5,width =0.1, color = 'k', alpha=0.8) #x,z,dx,dz
    ax.add_patch(arrowfig)
    art3d.pathpatch_2d_to_3d(arrowfig, z=0.6, zdir="y")
    
    if withoutbackground:
        ax.grid(False)
        ax.set(xlabel = '', ylabel='', zlabel='')
        ax.tick_params(left = False, right = False , labelleft = False , labelbottom = False, bottom = False)
        fig.patch.set_visible(False)
        ax.set_axis_off()
    
    if save:
        # Trim output
        bbox = fig.bbox_inches.from_bounds(2, 3.75, 12.25, 8) #AMOUNT TO CUT INTO FIGURE. left, top, width, height. Centered when left+0.5width=figurewidth, bottom+0.5height=figureheight
        plt.savefig('../output_images/bighero_with_arrows.pdf', bbox_inches=bbox)
    
    

#%%

plot(save=False, withoutbackground=False)


