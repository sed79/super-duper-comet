#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 15:04:21 2023

@author: sed79

Makes plots of S' comapring literature and MUSE data
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.legend_handler import HandlerTuple

perihelion_jd = 2458825.5
data = pd.read_csv("../extra/LiteratureS'forPlotting.csv")
dates = data['jd']-perihelion_jd
dataMUSE = pd.read_csv("../extra/MUSES'forPlotting.csv")
datesMUSE = dataMUSE['jd']-perihelion_jd

x = datesMUSE
y = dataMUSE["S"]
err = dataMUSE['uncertainty']

x_ = x[:-2]
y_ = y[:-2]
err_ = err[:-2]

print("UNWEIGHTED : EXCLUDING DATA AFTER OUTBURST ")
x_model = sm.add_constant(x_)
wls_model = sm.regression.linear_model.WLS(y_, x_model)#, weights=([1,1,1,1,1,5,1,1,1,1,1,1,1,1]))
results = wls_model.fit()
print(results.summary())
M = [results.params[1], results.params[0], f'Excluding data, unweighted. p-value: {results.pvalues[1]:.6f}']
predictions = results.get_prediction()
CI = predictions.conf_int(alpha=0.05)
CIpred = predictions.conf_int(alpha=0.05, obs=True)

CIcontinuous = []
CIxaxis = range(-30, 125)
for i in CIxaxis:   
    x0 = ([1, i])
    s = results.get_prediction(x0)
    interval = s.conf_int(alpha=0.05, obs=True)
    CIcontinuous.append([interval[0, 0], interval[0, 1]])
CIcontinuous = np.array(CIcontinuous)

def S_by_date_prediction():
    # Figure One
    fig, axes = plt.subplots(2, 1, figsize=(6,7), gridspec_kw={'height_ratios':[1, 0.6]})
    
    ax = axes[0]
    # Set color map and norm
    cmap = mpl.cm.Paired
    norm = mpl.colors.Normalize(vmin=400, vmax=1600)
    scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

    #Plotting
    p1 = ax.errorbar(dates, data["S"],yerr=data['uncertainty'], label='',  fmt='o', capsize=3,  c='grey', markersize=4, elinewidth=1)
    p2 = ax.scatter(dates, data["S"], zorder=5, cmap=cmap, norm=norm, c=data['wavelength midpoint'], s=80, linewidths=0.5, edgecolor='k')
    for i in range(len(data['jd'])):
        x = [(data['jd'][i] - data['date span'][i]-perihelion_jd), (data['jd'][i] + data['date span'][i]-perihelion_jd)]
        y = [data['S'][i], data['S'][i]]
        p3 = ax.plot(x, y, linestyle='-', linewidth=5, c=scalarMap.to_rgba(data['wavelength midpoint'][i]))
        
    p3 = ax.errorbar(datesMUSE, dataMUSE["S"],yerr=dataMUSE['uncertainty'], label='',  fmt='o', capsize=3,  c='grey', markersize=4, elinewidth=1)
    p4 = ax.scatter(datesMUSE, dataMUSE["S"], zorder=5, cmap=cmap, norm=norm, c=dataMUSE['wavelength midpoint'], s=80, marker='s', edgecolors='k', linewidths=0.5)
    
    p4 = ax.scatter([], [], c='w', edgecolor='k', linewidth=0.5, s=80, marker='o')
    p5 = ax.scatter([], [], c='w', edgecolor='k', linewidth=0.5, s=80, marker='s')
    
    fig.legend([p5, p4], ['MUSE','Other Work'], numpoints=1,  handler_map={tuple: HandlerTuple(ndivide=None)}, bbox_to_anchor=(0.8, 0.92))
    
    #Extras
    cb1 = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,fraction=0.1, shrink=0.8, location='right', pad=0.04, label='Wavelength range midpoint (nm)')    
    ax.set(xlabel='Days Relative to Perihelion', ylabel="S' (%/1000$\AA{}$)")
    ax.set_axisbelow(False)
    ax.grid(alpha=0.5, zorder=7)
    
    ax.set_ylim(0,27)
    
    # Top axis
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    new_tick_locations = [2458758,2458797,2458825,2458855,2458894,2458928] # JD corresponding to the Rh values in new_tick_labels
    new_tick_labels = ['-0.5', '-0.1', '0.0', '0.1','0.5','1.0']
    ax2.set_xticks(new_tick_locations - (np.ones_like(new_tick_locations)*perihelion_jd)+0.5)
    ax2.set_xticklabels(new_tick_labels)
    ax2.set_xlabel('Distance to Perihelion (au)')
    
    # ------------ second axis ---------------
    ax = axes[1]
    
    # shaded, dashed outline, white only 
    ax.fill_between(CIxaxis[0:116], CIcontinuous[0:116,1], CIcontinuous[0:116,0], color='tab:grey', alpha=0.4, linewidth=0, label='Prediction on new data points - model made without weights')
    ax.fill_between(CIxaxis[116:150], CIcontinuous[116:150,1]-0.1, CIcontinuous[116:150,0]+0.1, linewidth=2, color='tab:grey', linestyle='dashed', edgecolor='tab:grey', facecolor='white', alpha=0.5, label='Prediction on new data points - model made without weights')
    ax.fill_between(CIxaxis[114:], CIcontinuous[114:,1]-0.1, CIcontinuous[114:,0]+0.1, linewidth = 0, facecolor='white', alpha=1, zorder=2, label='Prediction on new data points - model made without weights')
    l0 = ax.vlines(x=300, ymin=-2, ymax=15, color='tab:grey', linestyle='dashed')
    #Plotting
    
    p3 = ax.errorbar(datesMUSE, dataMUSE["S"],yerr=dataMUSE['uncertainty'], label='',  fmt='o', capsize=3,  c='grey', markersize=4, elinewidth=1)
    p4 = ax.scatter(datesMUSE, dataMUSE["S"], zorder=5, color='white', norm=norm, s=80, marker='s', edgecolors='k', linewidths=0.5)
    
    p4 = ax.scatter([], [], c='w', edgecolor='k', linewidth=0.5, s=80, marker='o')
    p5 = ax.scatter([], [], c='w', edgecolor='k', linewidth=0.5, s=80, marker='s')
    
    #Extras
    ax.set(xlabel='Days Relative to Perihelion', ylabel="S' (%/1000$\AA{}$)")
    ax.set_axisbelow(False)
    ax.grid(alpha=0.5, zorder=7)
    
    ax.set_ylim(-1.5,15)
    ax.set_xlim(-35, 110)
    l1 = ax.vlines(x=2458913.50000-perihelion_jd, ymin=-2, ymax=15, color='xkcd:dark grey', linestyle='dashed')
    fig.legend([l0, l1], ['95% Predictive Interval', 'Outburst'], numpoints=1,  handler_map={tuple: HandlerTuple(ndivide=None)}, bbox_to_anchor=(0.5, 0.17))
    
    
    # Top axis
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    new_tick_locations = [2458797,2458825,2458855,2458867.5,2458886,2458901,2458915,2458928] # JD corresponding to the Rh values in new_tick_labels
    new_tick_labels = [ '-0.1', '0.0', '0.1', '0.2','0.4','0.6','0.8','1.0']

    ax2.set_xticks(new_tick_locations - np.array(perihelion_jd)+0.5)
    ax2.set_xticklabels(new_tick_labels)
    ax2.set_xlabel('Distance to Perihelion (au)')
    
    fig.tight_layout()
    plt.subplots_adjust(hspace=0.5)
    
    plt.savefig('../output_images/S_by_date_diff_5000km_prediction.pdf')
    return

def S_by_date():
    # Figure One
    fig, ax = plt.subplots(figsize=(6,4.5))
    # Set color map and norm
    cmap = mpl.cm.Paired
    norm = mpl.colors.Normalize(vmin=400, vmax=1600)
    scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    
    #Plotting
    p1 = ax.errorbar(dates, data["S"],yerr=data['uncertainty'], label='',  fmt='o', capsize=3,  c='grey', markersize=4, elinewidth=1)
    p2 = ax.scatter(dates, data["S"], zorder=5, cmap=cmap, norm=norm, c=data['wavelength midpoint'], s=80, linewidths=0.5, edgecolor='k')
    for i in range(len(data['jd'])):
        x = [(data['jd'][i] - data['date span'][i]-perihelion_jd), (data['jd'][i] + data['date span'][i]-perihelion_jd)]
        y = [data['S'][i], data['S'][i]]
        p3 = ax.plot(x, y, linestyle='-', linewidth=5, c=scalarMap.to_rgba(data['wavelength midpoint'][i]))
        
    p3 = ax.errorbar(datesMUSE, dataMUSE["S"],yerr=dataMUSE['uncertainty'], label='',  fmt='o', capsize=3,  c='grey', markersize=4, elinewidth=1)
    p4 = ax.scatter(datesMUSE, dataMUSE["S"], zorder=5, cmap=cmap, norm=norm, c=dataMUSE['wavelength midpoint'], s=80, marker='s', edgecolors='k', linewidths=0.5)
    
    p4 = ax.scatter([], [], c='w', edgecolor='k', linewidth=0.5, s=80, marker='o')
    p5 = ax.scatter([], [], c='w', edgecolor='k', linewidth=0.5, s=80, marker='s')
    
    fig.legend([p5, p4], ['MUSE','Other Work'], numpoints=1,  handler_map={tuple: HandlerTuple(ndivide=None)}, bbox_to_anchor=(0.78, 0.87))
    
    #Extras
    cb1 = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,fraction=0.1, shrink=0.8, location='right', pad=0.04, label='Wavelength range midpoint (nm)')    
    ax.set(xlabel='Days Relative to Perihelion', ylabel="S' (%/1000$\AA{}$)")
    ax.grid(alpha=0.5)
    
    # Top axis
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    new_tick_locations = [2458758,2458797,2458825,2458855,2458894,2458928] # JD corresponding to the Rh values in new_tick_labels
    new_tick_labels = ['-0.5', '-0.1', '0.0', '0.1','0.5','1.0']
    ax2.set_xticks(new_tick_locations - (np.ones_like(new_tick_locations)*perihelion_jd)+0.5)
    ax2.set_xticklabels(new_tick_labels)
    ax2.set_xlabel('Distance to Perihelion (au)')
    
    plt.savefig('../output_images/S_by_date.pdf')
    return

def S_by_wavelength_without_MUSE():
    fig, ax = plt.subplots(figsize=(6,8))
    # Set colormap and norm
    cmap = mpl.cm.RdBu_r
    norm = mpl.colors.CenteredNorm(vcenter=0)
    scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    # Plotting
    p1 = ax.errorbar(data['wavelength midpoint'], data["S"],yerr=data['uncertainty'], label='', 
                fmt='o', capsize=3,  c='grey', markersize=4)
    p2 = ax.scatter(data['wavelength midpoint'], data["S"], zorder=5, cmap=cmap, norm=norm, c=dates, edgecolor='dimgrey', linewidth=0.5)
    for i in range(len(data['jd'])):
        x = [(data['wavelength midpoint'][i] - data['wavelength range'][i]), (data['wavelength midpoint'][i] + data['wavelength range'][i])]
        y = [data['S'][i], data['S'][i]]
        p3 = ax.plot(x, y, linestyle='--', linewidth=1, c='grey')#scalarMap.to_rgba(dates[i]))
    
    # Extras
    cb1 = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,fraction=0.1, shrink=0.8, location='right', pad=0.04, label='Days Relative to Perihelion')    
    #ax.vlines(x=530, ymin=3, ymax=27, linestyle='dashed', color='k')
    #ax.vlines(x=860, ymin=3, ymax=27, linestyle='dashed', color='k')
    ax.set(xlabel='Wavelength (nm)', ylabel="S' (%/1000$\AA{}$)", ylim=(3,27))
    return

def S_by_wavelength():
    # Figure Three
    fig, ax = plt.subplots(figsize=(4,6))
    
    # Set colormap and norm
    cmap = mpl.cm.PuOr_r
    norm = mpl.colors.CenteredNorm(vcenter=0)
    scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    # Plotting Literature Values
    p1 = ax.errorbar(data['wavelength midpoint'], data["S"],yerr=data['uncertainty'], label='', 
                fmt='o', capsize=3,  c='grey', markersize=4, elinewidth=1)
    p2 = ax.scatter(data['wavelength midpoint'], data["S"], zorder=5, cmap=cmap, norm=norm, c=dates, edgecolor='k', linewidth=0.5, s=80, label='Other Work')
    for i in range(len(data['jd'])):
        x = [(data['wavelength midpoint'][i] - data['wavelength range'][i]), (data['wavelength midpoint'][i] + data['wavelength range'][i])]
        y = [data['S'][i], data['S'][i]]
        p3 = ax.plot(x, y, linestyle='--', linewidth=1, c='grey')#scalarMap.to_rgba(dates[i]))
    # Plotting MUSE Values
    p1 = ax.errorbar(dataMUSE['wavelength midpoint'], dataMUSE["S"],yerr=dataMUSE['uncertainty'], label='', 
                fmt='o', capsize=3,  c='grey', markersize=4, elinewidth=1)
    p2 = ax.scatter(dataMUSE['wavelength midpoint'], dataMUSE["S"], zorder=5, cmap=cmap, norm=norm, c=datesMUSE, edgecolor='k', linewidth=0.5, s=80, marker='s', label='MUSE')
    for i in range(len(dataMUSE['jd'])):
        x = [(dataMUSE['wavelength midpoint'][i] - dataMUSE['wavelength range'][i]), (dataMUSE['wavelength midpoint'][i] + dataMUSE['wavelength range'][i])]
        y = [dataMUSE['S'][i], dataMUSE['S'][i]]
        p3 = ax.plot(x, y, linestyle='--', linewidth=1, c='grey')#c=scalarMap.to_rgba(datesMUSE[i]))
    
    p4 = ax.scatter([], [], c='w', edgecolor='k', linewidth=0.5, s=80, marker='o')
    p5 = ax.scatter([], [], c='w', edgecolor='k', linewidth=0.5, s=80, marker='s')
    
    # Extras
    cb1 = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,fraction=0.1, shrink=0.8, location='right', pad=0.04, label='Days Relative to Perihelion')    
    cb1.ax.invert_yaxis()
    
    ax.set(xlabel='Wavelength (nm)', ylabel="S' (%/1000$\AA{}$)", ylim=(3,27))
    ax.grid(alpha=0.5)
    fig.legend([p5, p4], ['MUSE','Other Work'], numpoints=1,  handler_map={tuple: HandlerTuple(ndivide=None)}, bbox_to_anchor=(0.75, 0.97))
    plt.tight_layout()
    plt.savefig('../output_images/S_by_lambda.pdf')
    return

#S_by_date()
#S_by_wavelength()
S_by_date_prediction
