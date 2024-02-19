#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 16:03:51 2023

@author: sed79

Plots the production rates of 2I using literature values from 
../extra/otherprodrates.csv and MUSE data from ../extra/MUSEprodrates.csv
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
from matplotlib.transforms import Bbox
from matplotlib import rc

import statsmodels.api as sm
import statsmodels


def plot():

    #GENERATE CONFIDENCE INTERVALS
    x = data['jd'][4:]
    x_ = data['jd'][4:-2]
    CIaxis = np.linspace(x.iloc[0]-4, x.iloc[-1]+4, 50) # extend each side of the data range for overlap with marker
    CoverCy = data['ratio Q(C2)/Q(CN)'][4:]
    CoverCy_ = data['ratio Q(C2)/Q(CN)'][4:-2]
    NoverCy = data['ratio Q(NH2)/Q(CN)'][4:]
    NoverCy_ = data['ratio Q(NH2)/Q(CN)'][4:-2]

    CoverCyerr = data['ratio Q(C2)/Q(CN) un'][4:]
    CoverCyerr_ = data['ratio Q(C2)/Q(CN) un'][4:-2]
    NoverCyerr = data['ratio Q(NH2)/Q(CN) un'][4:]
    NoverCyerr_ = data['ratio Q(NH2)/Q(CN) un'][4:-2]
    
    y = CoverCy
    y_ = CoverCy_
    yerr = CoverCyerr
    yerr_ = CoverCyerr_

    x_model = sm.add_constant(x_)
    wls_model = sm.regression.linear_model.WLS(y_, x_model, weights=[1/x**2 for x in yerr_]) #err_
    results = wls_model.fit()
    print(results.summary())
    MCoverC_ = [results.params[1], results.params[0], f'C2/CN. p-value: {results.pvalues[1]:.6f}']
    # predictions = results.get_prediction()
    # CI4 = predictions.conf_int(alpha=0.05)
    # CI4pred = predictions.conf_int(alpha=0.05, obs=True)

    CICoverCcontinuous = []
    for i in CIaxis:   
        x0 = ([1, i])
        s = results.get_prediction(x0)
        interval = s.conf_int(alpha=0.05, obs=False) #obs True/False returns predictive/mean confidence interval. Changes betwee using standard error of observations or of mean.
        CICoverCcontinuous.append([interval[0, 0], interval[0, 1]])
    CICoverCcontinuous = np.array(CICoverCcontinuous)

    y = NoverCy
    y_ = NoverCy_
    yerr = NoverCyerr
    yerr_ = NoverCyerr_

    x_model = sm.add_constant(x_)
    wls_model = sm.regression.linear_model.WLS(y_, x_model, weights=[1/x**2 for x in yerr_]) #, weights=1/(yerr_)) #err_
    results = wls_model.fit()
    print(results.summary())
    MNoverC_ = [results.params[1], results.params[0], f'NH2/CN. p-value: {results.pvalues[1]:.6f}']

    CINoverCcontinuous = []
    for i in CIaxis:   
        x0 = ([1, i])
        s = results.get_prediction(x0)
        interval = s.conf_int(alpha=0.05, obs=False) #obs True/False returns predictive/mean confidence interval. Changes betwee using standard error of observations or of mean.
        CINoverCcontinuous.append([interval[0, 0], interval[0, 1]])
    CINoverCcontinuous = np.array(CINoverCcontinuous)
    
    # MAKE FIGURE

    data = pd.read_csv("../extra/MUSEprodrates.csv")
    otherdata = pd.read_csv("../extra/otherprodrates.csv")
    
        fig, axs = plt.subplots(3, 1, figsize=(14,16),gridspec_kw={'hspace': 0.04, 'height_ratios':[1.4, 1, 1]}, sharex=True)
    ax = axs[0]

    # -------------
    # FIRST SUBPLOT
    # -------------
    
    # C2 data
    p1 = ax.errorbar(data['jd'], data['Q(C2)']*10**24, yerr=data['Q(C2) un']*10**24, label='C$_2$ - this work', 
                fmt='s', capsize=3, markerfacecolor='#006755', markeredgecolor='#002D25', c='grey', markersize=10, elinewidth=1) #006755
    p2 = ax.errorbar(otherdata['jd'], otherdata['Q(C2)']*10**24, yerr=otherdata['Q(C2) un']*10**24, label='C$_2$ - other work', 
                fmt='o', capsize=3, markerfacecolor='w', markeredgecolor='#006755', c='grey', markersize=10, elinewidth=1)
    p3 = ax.errorbar(otherdata['jd'], otherdata['Q(C2) upper lim']*10**24, 
                fmt='v', markerfacecolor='w', markeredgecolor='#006755',  markersize=11, elinewidth=1)
    
    # Confidence interval - C2
    M = M5
    CI = CI5continuous
    # ax.plot(CIaxis, CI[:,1],  color='#002D25', linestyle='dashed', alpha=0.5, label='Predictive confidence interval') # top line
    # ax.plot(CIaxis, CI[:,0],  color='#002D25', linestyle='dashed', alpha=0.5, label='Predictive confidence interval') # bottom line

    
    # CN data
    p4 = ax.errorbar(data['jd'], data['Q(CN)']*10**24, yerr=data['Q(CN) un']*10**24, label='CN - this work', 
                fmt='s', capsize=3, markerfacecolor='#1E88E5', markeredgecolor='#0C3B65', c='grey', markersize=10, elinewidth=1, zorder=5) #1E88E5
    p5 = ax.errorbar(otherdata['jd'], otherdata['Q(CN)']*10**24, yerr=otherdata['Q(CN) un']*10**24, label='CN - other work', 
                fmt='o', capsize=3, markerfacecolor='w', markeredgecolor='#1E88E5', c='grey', markersize=10, elinewidth=1)
    p6 = ax.errorbar(otherdata['jd'], otherdata['Q(CN) upper lim']*10**24, 
                fmt='v', markerfacecolor='w', markeredgecolor='#1E88E5', markersize=11, elinewidth=1) #, yerr=c2_lower_error
    
    # Confidence interval - CN
    M = M6
    CI = CI6continuous
    # ax.fill_between(CIaxis ,[x*10**24 for x in CI[:,1]], [x*10**24 for x in CI[:,0]], color='tab:blue', alpha=0.5, label='Predictive confidence interval')
    # ax.plot(CIaxis ,CI[:,1],  color='#1E88E5', linestyle='dashed', alpha=0.5, label='Predictive confidence interval') # top line
    # ax.plot(CIaxis ,CI[:,0],  color='#1E88E5', linestyle='dashed', alpha=0.5, zorder=12, label='Predictive confidence interval') # bottom line

    # ax.fill_between(CIaxis[1:39], CI[1:39,1], CI[1:39,0], color='#1E88E5', alpha=0.3, linewidth=0, zorder=4, label='Prediction on new data points - model made without weights')
    # ax.fill_between(CIaxis[38:-1], CI[38:-1,1], CI[38:-1,0], linewidth=2, color='#1E88E5', linestyle='dashed', edgecolor='tab:grey', zorder=1, facecolor='white', alpha=0.5, label='') # reduced on x axis to have edges covered
    # ax.fill_between(CIaxis[25:], CI[25:,1]-0.02, CI[25:,0]+0.02, linewidth = 0, facecolor='white', alpha=1, zorder=1, label='Prediction on new data points - model made without weights') # reduced on y-axis to show dashed lines

    # NH2 data
    p7 = ax.errorbar(data['jd'], data['Q(NH2)']*10**24, yerr=data['Q(NH2) un']*10**24, label='NH$_2$ - this work', 
                fmt='s', capsize=3,  markerfacecolor='#D81B60', markeredgecolor='#550925', c='grey', markersize=10, elinewidth=1) 
    
    # Confidence interval - NH2
    M = M4
    CI = CI4continuous
    # #l = ax.plot(xx, (M[0]*xx + M[1])*10**24, label=M[2], color='tab:blue')
    # #ax.fill_between(np.linspace(x.iloc[0], x.iloc[-1], 20),[x*10**24 for x in CI[:,1]], [x*10**24 for x in CI[:,0]], color='tab:pink', alpha=0.5, label='Predictive confidence interval')
    # ax.plot(CIaxis ,CI[:,1],  color='#550925', linestyle='dashed', alpha=0.5, label='Predictive confidence interval') # top line
    # ax.plot(CIaxis ,CI[:,0],  color='#550925', linestyle='dashed', alpha=0.5, label='Predictive confidence interval') # bottom line

    # Perihelion and outburst
    p8 = ax.vlines(x=2458825, ymin=0.8e22, ymax=1e25, linestyle='--', linewidth=1, color='k', label='Perihelion')
    p9 = ax.vlines(x=2458915, ymin=0.8e22, ymax=1e25, linestyle='dotted', linewidth=1, color='k', label='Splitting')
    
    p10 = ax.scatter([], [], c='w', label="Extra label on the legend")
    
    font = {'family': 'DejaVu Sans',
        'weight': 'normal',
        'size': 10,
        }
    
    l = fig.legend([p1, p4, p7], 
                   ['C$_2$                   ', 'CN', 'NH$_2$'], #set number of spaces in label strings are important for alignment!
                   title='MUSE:      ', title_fontproperties=font, 
                    handler_map={tuple: HandlerTuple(ndivide=None)},
                   bbox_to_anchor=(1.03, 0.87), framealpha = 0) #figure (x,y) position, width, height
    
    l = fig.legend([p2, p5], 
                   ['C$_2$                   ', 'CN', 'NH$_2$'],
                   title='Other Work:', title_fontproperties=font, 
                  handler_map={tuple: HandlerTuple(ndivide=None)},
                  bbox_to_anchor=(1.03, 0.805), framealpha = 0) #figure (x,y) position
    
    l = fig.legend([p10], 
                   ['                        '], numpoints=1, 
                   title='Other Work', title_fontproperties=font,
                  handler_map={tuple: HandlerTuple(ndivide=None)},
                  bbox_to_anchor=(1.03, 0.75), framealpha = 0
                  ) 
    
    l = fig.legend([p3, p6], 
                   ['C$_2$                   ', 'CN'], numpoints=1,
                   title='   Upper Limits:', title_fontproperties=font,
                  handler_map={tuple: HandlerTuple(ndivide=None)},
                  bbox_to_anchor=(1.03, 0.738), framealpha = 0) #figure (x,y) position
    
    l = fig.legend([p8, p9], 
                   ['Perihelion       ', 'Outburst'], numpoints=1,
                  handler_map={tuple: HandlerTuple(ndivide=None)},
                  bbox_to_anchor=(1.03, 0.65), framealpha = 0) #figure (x,y) position

    
    xlocs = [2458757, 2458788, 2458818, 2458849, 2458880, 2458909, 2458940]
    xlabs = ['1 Oct 2019', '1 Nov 2019', '1 Dec 2019', '1 Jan 2020', '1 Feb 2020', '1 Mar 2020', '1 Apr 2020']
    
    ax.grid(alpha=0.5)
    ax.grid(which='minor', axis='y', alpha=0.5)
    ax.set(ylim=(1e23, 1e25), xticks=[])
    ax.set_ylabel('Production Rate (mol s$^{-1}$)', fontsize=14, labelpad=6)
    ax.tick_params(axis='x', color='white')
    
    # Top axis
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    new_tick_locations = [2458737,2458744,2458751,2458758,2458766,2458775,2458785,2458797,2458806,2458825,2458846,2458855,2458867.5,2458877,2458886,2458894,2458901,2458908,2458915,2458921,2458928,2458934,2458940] # JD corresponding to the Rh values in new_tick_labels
    new_tick_labels = ['-0.8', '-0.7', '-0.6', '-0.5', '-0.4', '-0.3', '-0.2', '-0.1', '-0.05', '0.0', '0.05', '0.1', '0.2', '0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','1.1','1.2']
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(new_tick_labels)
    ax2.set_xlabel('Distance to Perihelion (au)', fontsize=14, labelpad=10)
    
    ax.set_yscale('log')

    # -------------
    # SECOND SUBPLOT
    # -------------
    
    # --------------------------------------
    
    ax2 = axs[1]
    # Our data - C2 is 10^23 and CN is 10^24, hence x10 for CN data
    p1 = ax2.errorbar(data['jd'], data['ratio Q(C2)/Q(CN)'], yerr=data['ratio Q(C2)/Q(CN) un'], markerfacecolor='dimgrey', 
                markeredgecolor='k', c='grey', markersize=10, fmt='s', label='This work', capsize=3, elinewidth=1)
    
    # Other data - with both error bars
    p2 = ax2.errorbar(otherdata['jd'], otherdata['ratio Q(C2)/Q(CN)'], yerr=otherdata['ratio Q(C2)/Q(CN) un'], markerfacecolor='w', 
                markeredgecolor='k', c='grey', markersize=10, fmt='o', label='Other work', capsize=3, elinewidth=1)
    
    # Other data - With upper error shown by v symbol
    p3 = ax2.errorbar(otherdata['jd'], otherdata['ratio Q(C2)/Q(CN) upper lim'], markerfacecolor='w', 
                markeredgecolor='k', c='grey', markersize=11, fmt='v', label='Other work', capsize=3, elinewidth=1)
    
    # Confidence Interval
    M = MCoverC_
    CI = CICoverCcontinuous
    #ax.fill_between(CIaxis,CI[:,1],CI[:,0], color='tab:grey', alpha=0.5, label='Predictive confidence interval')
    p11 = ax2.plot(CIaxis ,CI[:,1],  color='tab:grey', linestyle='dashed', alpha=0.5, label='Predictive confidence interval') # top line
    ax2.plot(CIaxis ,CI[:,0],  color='tab:grey', linestyle='dashed', alpha=0.5, label='Predictive confidence interval') # bottom line

    # Perihelion
    p8 = ax2.vlines(x=2458825, ymin=0, ymax=0.8, linestyle='--', linewidth=1, color='k', label='Perihelion')
    p9 = ax2.vlines(x=2458915, ymin=0, ymax=0.8, linestyle='dotted', linewidth=1, color='k', label='Splitting')
    
    # Dates x axis
    ax2.tick_params(axis='x', color='white')
    ax2.set(xlabel='', 
           ylim=(0, 0.8), xlim=(xlocs[0]-20, xlocs[-1]))
    ax2.set_ylabel('Q(C$_2$)/Q(CN)', fontsize=14, labelpad=8)
    
    # Grid and legend
    ax2.grid(alpha=0.5)
    l = fig.legend([p1, p2, p3], 
                   ['MUSE              ', 'Other Work', 'Other Work'], numpoints=1, #set number of spaces in label strings are important for alignment!
                  handler_map={tuple: HandlerTuple(ndivide=None)},
                  bbox_to_anchor=(1.03, 0.54), framealpha = 0) #(1.02, 0.3) for skinny legend
    
    l = fig.legend([p10], 
                   ['        Upper Limit     '], numpoints=1, 
                  handler_map={tuple: HandlerTuple(ndivide=None)},
                  bbox_to_anchor=(1.03, 0.503), framealpha = 0) 
    
    l = fig.legend([p8, p9], 
                   ['Perihelion       ', 'Outburst'], numpoints=1,
                  handler_map={tuple: HandlerTuple(ndivide=None)},
                  bbox_to_anchor=(1.03, 0.47), framealpha = 0) #figure (x,y) position
    
    l = fig.legend([p11[0]], 
                ['Q(C$_2$)/Q(CN)'], numpoints=1,
                title=' 95% prediction \nconfidence interval:', title_fontproperties=font,
                handler_map={tuple: HandlerTuple(ndivide=None)},
                bbox_to_anchor=(1.03, 0.43), framealpha = 0) #figure (x,y) position
    
    # Remove righthand side edge
    #axs[1].spines[['right']].set_visible(False)
    
    # -------------
    # THIRD SUBPLOT
    # -------------
    
    #-------------------------------------------------------------------------
    
    ax3 = axs[2]
    # Our data - both NH2 and CN in units 10^24 mol/sec
    p1 = ax3.errorbar(data['jd'], data['ratio Q(NH2)/Q(CN)'], yerr=data['ratio Q(NH2)/Q(CN) un'], markerfacecolor='dimgrey', 
                markeredgecolor='k', c='grey', markersize=10, fmt='s', label='This work', capsize=3, elinewidth=1, zorder=5)
    
    # Confidence Interval
    M = MNoverC_
    CI = CINoverCcontinuous
    #p11 = ax3.fill_between(CIaxis[1:39], CI[1:39,1], CI[1:39,0], color='tab:grey', alpha=0.4, linewidth=0, zorder=4, label='Prediction on new data points - model made without weights')
    p11 = ax3.plot(CIaxis ,CI[:,1],  color='tab:grey', linestyle='dashed', alpha=0.5, label='Predictive confidence interval') # top line
    ax3.plot(CIaxis ,CI[:,0],  color='tab:grey', linestyle='dashed', alpha=0.5, label='Predictive confidence interval') # bottom line

    # Perihelion
    p8 = ax3.vlines(x=2458825, ymin=0, ymax=3, linestyle='--', linewidth=1, color='k', label='Perihelion')
    p9 = ax3.vlines(x=2458915, ymin=0, ymax=3, linestyle='dotted', linewidth=1, color='k', label='Splitting')
    
    # Dates x axis
    xlocs = [2458757, 2458788, 2458818, 2458849, 2458880, 2458909, 2458940]
    xlabs = ['1 Oct 2019', '1 Nov 2019', '1 Dec 2019', '1 Jan 2020', '1 Feb 2020', '1 Mar 2020', '1 Apr 2020']
    ax3.set(xticks=xlocs, xticklabels=xlabs, 
           xlabel='', 
           ylim=(0, 3), xlim=(xlocs[0]-20, xlocs[-1]))
    ax3.set_ylabel('Q(NH$_2$)/Q(CN)', fontsize=14, labelpad=8)
    
    # Grid and legend
    ax3.grid(alpha=0.5, zorder=10)
    l = fig.legend([p1], 
                   ['MUSE              '], numpoints=1, #set number of spaces in label strings are important for alignment!
                  handler_map={tuple: HandlerTuple(ndivide=None)},
                  bbox_to_anchor=(1.03, 0.28), framealpha = 0) #(1.02, 0.3) for skinny legend
    
    
    l = fig.legend([p8, p9], 
                   ['Perihelion       ', 'Outburst'], numpoints=1,
                  handler_map={tuple: HandlerTuple(ndivide=None)},
                  bbox_to_anchor=(1.03, 0.25), framealpha = 0) #figure (x,y) position
    
    l = fig.legend([p11[0]], 
                ['Q(NH$_2$)/Q(CN)'], numpoints=1,
                title=' 95% prediction \nconfidence interval:', title_fontproperties=font,
                handler_map={tuple: HandlerTuple(ndivide=None)},
                bbox_to_anchor=(1.03, 0.21), framealpha = 0) #figure (x,y) position
    
    bbox = (fig.bbox_inches.from_bounds(1, 1.5, 13.4, 13.25)) #perfect numbers are (1, 1, 13.5, 10.25)
    plt.savefig('../figures/ProductionRates&Ratios_prediction.pdf',bbox_inches=bbox)
    
    
    return

plot()


