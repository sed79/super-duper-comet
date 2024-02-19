#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
from matplotlib.transforms import Bbox
from matplotlib import rc

import statsmodels.api as sm
import statsmodels

perihelion_jd = 2458825.5
dates = data['jd']-perihelion_jd

dataMUSE = pd.read_csv("../extra/MUSES'forPlotting_diff_5000km.csv")
datesMUSE = dataMUSE['jd']-perihelion_jd

data = pd.read_csv("../extra/MUSEprodrates.csv")
otherdata = pd.read_csv("../extra/otherprodrates.csv")

x = datesMUSE
y = dataMUSE["S"]
err = dataMUSE['uncertainty']

# data pre outburst
x_ = x[:-2]
y_ = y[:-2]
err_ = err[:-2]

# fit to data excluding post outburst
x_model = sm.add_constant(x_)
wls_model = sm.regression.linear_model.WLS(y_, x_model, weights=([1/x**2 for x in err_]))
results = wls_model.fit()
print(results.summary())
M = [results.params[1], results.params[0], f'Excluding data, unweighted. p-value: {results.pvalues[1]:.6f}']
predictions = results.get_prediction()
# CI = predictions.conf_int(alpha=0.05)
# CIpred = predictions.conf_int(alpha=0.05, obs=True)

CIcontinuous = []
CIxaxisS = np.array(range(-30, 115))
for i in CIxaxisS:   
    x0 = ([1, i])
    s = results.get_prediction(x0)
    interval = s.conf_int(alpha=0.05, obs=False)
    CIcontinuous.append([interval[0, 0], interval[0, 1]])
CIcontinuous = np.array(CIcontinuous)

gridspec = dict(hspace=0.0, height_ratios=[1, 0.04, 1, 0.04, 1, 0.5, 1])
fig, axs = plt.subplots(nrows=7, ncols=1, figsize=(5,15), gridspec_kw=gridspec)

ax1, ax0, ax2, ax00, ax3, ax000, ax4 = axs

axs = [ax1, ax2, ax3, ax4]

ax0.set_visible(False)
ax00.set_visible(False)
ax000.set_visible(False)

# C2 data
p1 = ax1.errorbar(data['jd'][4:], data['Q(C2)'][4:]*10**24, yerr=data['Q(C2) un'][4:]*10**24, label='C$_2$ - this work', # all but first 2 data points
            fmt='s', capsize=3, markerfacecolor='#006755', markeredgecolor='#002D25', c='grey', markersize=10, elinewidth=1) #006755


# Confidence interval - C2
M = M5
CI = CI5continuous

ax1.plot(CIaxis ,CI[:,1],  color='#006755', linestyle='dashed', alpha=0.5, label='Predictive confidence interval') # top line
p11 = ax1.plot(CIaxis ,CI[:,0],  color='#006755', linestyle='dashed', alpha=0.5, label='Predictive confidence interval') # bottom line
ax1.fill_between(CIaxis[:39], CI[:39,1], CI[:39,0], color='#006755', alpha=0.3, linewidth=0, zorder=4, label='Prediction on new data points - model made with weights')


#ax1.fill_between(CIaxis[:39], CI[:39,1], CI[:39,0], color='#006755', alpha=0.4, linewidth=0, zorder=4, label='Prediction on new data points - model made without weights')
#ax1.fill_between(CIaxis[38:-1], CI[38:-1,1], CI[38:-1,0], linewidth=3, color='#002D25', linestyle='dashed', edgecolor='tab:grey', zorder=1, facecolor='white', alpha=0.8, label='') # reduced on x axis to have edges covered
#ax1.fill_between(CIaxis[25:], CI[25:,1]-0.02, CI[25:,0]+0.02, linewidth = 0, facecolor='white', alpha=1, zorder=1, label='Prediction on new data points - model made without weights') # reduced on y-axis to show dashed lines


# CN data
p4 = ax3.errorbar(data['jd'][4:], data['Q(CN)'][4:]*10**24, yerr=data['Q(CN) un'][4:]*10**24, label='CN - this work', 
            fmt='s', capsize=3, markerfacecolor='#1E88E5', markeredgecolor='#0C3B65', c='grey', markersize=10, elinewidth=1, zorder=5) #1E88E5

# Confidence interval - CN
M = M6
CI = CI6continuous
#ax3.fill_between(CIaxis ,[x*10**24 for x in CI[:,1]], [x*10**24 for x in CI[:,0]], color='tab:blue', alpha=0.5, label='Predictive confidence interval')
p12 = ax3.plot(CIaxis ,CI[:,1],  color='#1E88E5', linestyle='dashed', alpha=0.8, label='Predictive confidence interval') # top line
ax3.plot(CIaxis ,CI[:,0],  color='#1E88E5', linestyle='dashed', alpha=0.8, zorder=12, label='Predictive confidence interval') # bottom line
ax3.fill_between(CIaxis[:39], CI[:39,1], CI[:39,0], color='#1E88E5', alpha=0.3, linewidth=0, zorder=4, label='Prediction on new data points - model made with weights')

# ax3.fill_between(CIaxis[1:39], CI[1:39,1], CI[1:39,0], color='#1E88E5', alpha=0.3, linewidth=0, zorder=4, label='Prediction on new data points - model made without weights')
# ax3.fill_between(CIaxis[38:-1], CI[38:-1,1], CI[38:-1,0], linewidth=2, color='#1E88E5', linestyle='dashed', edgecolor='tab:grey', zorder=1, facecolor='white', alpha=0.5, label='') # reduced on x axis to have edges covered
# ax3.fill_between(CIaxis[25:], CI[25:,1]-0.02, CI[25:,0]+0.02, linewidth = 0, facecolor='white', alpha=1, zorder=1, label='Prediction on new data points - model made without weights') # reduced on y-axis to show dashed lines

# NH2 data
p7 = ax2.errorbar(data['jd'][4:], data['Q(NH2)'][4:]*10**24, yerr=data['Q(NH2) un'][4:]*10**24, label='NH$_2$ - this work', 
            fmt='s', capsize=3,  markerfacecolor='#D81B60', markeredgecolor='#550925', c='grey', markersize=10, elinewidth=1) 

# Confidence interval - NH2
M = M4
CI = CI4continuous
#l = ax.plot(xx, (M[0]*xx + M[1])*10**24, label=M[2], color='tab:blue')
#ax.fill_between(np.linspace(x.iloc[0], x.iloc[-1], 20),[x*10**24 for x in CI[:,1]], [x*10**24 for x in CI[:,0]], color='tab:pink', alpha=0.5, label='Predictive confidence interval')
ax2.plot(CIaxis ,CI[:,1],  color='#550925', linestyle='dashed', alpha=0.5, label='Predictive confidence interval') # top line
p13 = ax2.plot(CIaxis ,CI[:,0],  color='#550925', linestyle='dashed', alpha=0.5, label='Predictive confidence interval') # bottom line
ax2.fill_between(CIaxis[:39], CI[:39,1], CI[:39,0], color='#D81B60', alpha=0.3, linewidth=0, zorder=4, label='Prediction on new data points - model made with weights')
ax2.fill_between(CIaxis[38:-1], CI[38:-1,1], CI[38:-1,0], linewidth=2, color='#D81B60', linestyle='dashed', edgecolor='tab:grey', zorder=1, facecolor='white', alpha=0.5, label='') # reduced on x axis to have edges covered
ax2.fill_between(CIaxis[25:], CI[25:,1]-0.02, CI[25:,0]+0.02, linewidth = 0, facecolor='white', alpha=1, zorder=1, label='Prediction on new data points - model made with weights') # reduced on y-axis to show dashed lines

# DUST COLOUR

# shading and dashed lines
ax4.plot(CIxaxisS+perihelion_jd ,CIcontinuous[:,1],  color='tab:grey', linestyle='dashed', alpha=0.5, label='Predictive confidence interval') # top line
p15 = ax4.plot(CIxaxisS+perihelion_jd ,CIcontinuous[:,0],  color='tab:grey', linestyle='dashed', alpha=0.5, label='Predictive confidence interval') # bottom line
ax4.fill_between(CIxaxisS[0:116]+perihelion_jd, CIcontinuous[0:116,1], CIcontinuous[0:116,0], color='tab:grey', alpha=0.4, linewidth=0, label='Prediction on new data points - model made without weights')
#ax4.fill_between(CIxaxisS[116:150]+perihelion_jd, CIcontinuous[116:150,1]-0.1, CIcontinuous[116:150,0]+0.1, linewidth=2, color='tab:grey', linestyle='dashed', edgecolor='tab:grey', facecolor='white', alpha=0.5, label='Prediction on new data points - model made without weights')
#ax4.fill_between(CIxaxisS[114:]+perihelion_jd, CIcontinuous[114:,1]-0.1, CIcontinuous[114:,0]+0.1, linewidth = 0, facecolor='white', alpha=1, zorder=2, label='Prediction on new data points - model made without weights')
l0 = ax4.vlines(x=300, ymin=-2, ymax=15, color='tab:grey', linestyle='dashed')

# data
p3 = ax4.errorbar(datesMUSE+perihelion_jd, dataMUSE["S"],yerr=dataMUSE['uncertainty'], label='',  fmt='o', capsize=3,  c='grey', markersize=4, elinewidth=1)
p4 = ax4.scatter(datesMUSE+perihelion_jd, dataMUSE["S"], zorder=5, color='white', s=80, marker='s', edgecolors='k', linewidths=0.5)

p4 = ax4.scatter([], [], c='w', edgecolor='k', linewidth=0.5, s=80, marker='o')
p5 = ax4.scatter([], [], c='w', edgecolor='k', linewidth=0.5, s=80, marker='s')


#Extras
ax.set(xlabel='Days Relative to Perihelion', ylabel="S' (%/1000$\AA{}$)")
ax.set_axisbelow(False)
ax.grid(alpha=0.5, zorder=7)

# ax.set_ylim(-1.5,15)
# ax.set_xlim(-35, 110)
# l1 = ax.vlines(x=2458915 - perihelion_jd, ymin=-2, ymax=15, color='k', linestyle='dotted', linewidth=1)

# Outburst
for ax in axs:
    p9 = ax.vlines(x=2458915, ymin=0.8e22, ymax=1e25, linestyle='dotted', linewidth=1, color='k', label='Splitting')
ax4.vlines(x=2458915, ymin=-5, ymax=25, linestyle='dotted', linewidth=1, color='k', label='Splitting')

# Perihelion
#p8 = ax4.vlines(x=2458825, ymin=ax4.get_ylim()[0], ymax=ax4.get_ylim()[1], linestyle='--', linewidth=1, color='k', label='Perihelion')

p10 = ax1.plot([1,2],[1,2],  color='k', linestyle='dashed', alpha=0.5, label='created for legend')
#p14 = ax1.arrow(0,0,-30,0, fc='w', ec='k', alpha=0.9, head_width=0.2e24,head_length=10,width=0.1e24,length_includes_head=True)
p14 = ax.scatter(1,1, c='k', marker=r'$\leftarrow$',s=150)

font = {'family': 'DejaVu Sans',
    'weight': 'normal',
    'size': 10,
    }

# l = fig.legend([p1, p7, p4], 
#                 ['C$_2$                   ',  'NH$_2$', 'CN'], #set number of spaces in label strings are important for alignment!
#                 title='MUSE:      ', title_fontproperties=font, 
#                 handler_map={tuple: HandlerTuple(ndivide=None)},
#                 bbox_to_anchor=(1.03, 0.87), framealpha = 0) #figure (x,y) position, width, height
l = fig.legend([p14], 
                [' \nHST photometry\ndates'], numpoints=1,
                handler_map={tuple: HandlerTuple(ndivide=None)},
                bbox_to_anchor=(1.22, 0.895), framealpha = 0) #figure (x,y) position

for ypos, p_outburst in zip([0.14, 0.32, 0.5], [p12, p13, p11]):

    l = fig.legend([p9], 
                    ['Outburst           '], numpoints=1,
                    handler_map={tuple: HandlerTuple(ndivide=None)},
                    bbox_to_anchor=(1.26, ypos+0.34), framealpha = 0) #figure (x,y) position

    l = fig.legend([p_outburst[0]], 
                    ['Prediction \nConfidence \nInterval 95%     '],
                    handler_map={tuple: HandlerTuple(ndivide=None)},
                    bbox_to_anchor=(1.26, ypos+0.32), framealpha = 0) #figure (x,y) position

# l = fig.legend([p8], 
#                     ['Perihelion         '], numpoints=1,
#                     handler_map={tuple: HandlerTuple(ndivide=None)},
#                     bbox_to_anchor=(1.26, 0.23), framealpha = 0) #figure (x,y) position
    
l = fig.legend([p9], 
                    ['Outburst           '], numpoints=1,
                    handler_map={tuple: HandlerTuple(ndivide=None)},
                    bbox_to_anchor=(1.26, 0.22), framealpha = 0) #figure (x,y) position

l = fig.legend([p10[0]], 
                ['Prediction \nConfidence \nInterval 95%     '],
                handler_map={tuple: HandlerTuple(ndivide=None)},
                bbox_to_anchor=(1.26, 0.20), framealpha = 0) #figure (x,y) position


# Set y axis limits and labels
ax1.set(ylim=(0.1e22, 0.9e24), yticks=[0.1e24, 0.2e24, 0.3e24, 0.4e24, 0.5e24, 0.6e24, 0.7e24, 0.8e24, 0.9e24],
        yticklabels=['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'])
ax2.set(ylim=(1e23, 0.7e25), yticks=[1e24, 2e24, 3e24, 4e24, 5e24, 6e24, 7e24],
        yticklabels=['1', '2', '3', '4', '5', '6', '7'])
ax3.set(ylim=(0.5e24, 0.4e25), yticks=[0.5e24, 1e24,1.5e24, 2e24, 2.5e24, 3e24, 3.5e24, 4e24],
        yticklabels=['0.5', '1.0', '1.5', '2.0', '2.5', '3.0', '3.5', '4.0'])
ax4.set(ylim=(0,15))




# Inner ticks
for jewitt_date in [2458851.5, 2458877.5, 2458903.5, 2458931.5, 2458938.5, 2458945.5, 2458952.5, 2458959.5]:
        ax1.arrow(jewitt_date, 0.9e24, 0, -0.03e24, fc='k', ec='k', alpha=0.9, head_width=0,head_length=0,width=0.01,length_includes_head=True)
        # ax2.arrow(jewitt_date, 7e24, 0, -0.25e24, fc='k', ec='k', alpha=0.9, head_width=0,head_length=0,width=0.01,length_includes_head=True)
        # ax3.arrow(jewitt_date, 4.0e24, 0, -0.15e24, fc='k', ec='k', alpha=0.9, head_width=0,head_length=0,width=0.01,length_includes_head=True)

# Labels
props = dict(boxstyle='square', facecolor='white', alpha=0.7, pad=0.2, edgecolor='k')
ax1.text(2458830.5, 0.06e24,'C$_2$', color='k',fontsize=18,bbox=props)#, fontweight='bold')
ax2.text(2458830.5, 0.55e24,'NH$_2$', color='k',fontsize=18,bbox=props)
ax3.text(2458829.5, 0.7e24,'CN', color='k',fontsize=18,bbox=props)
ax4.text(2458831-40, 1,"S'", color='k',fontsize=18,bbox=props)

for j, ax in enumerate(axs):
    ax.grid(alpha=0.5)
    ax.grid(which='minor', axis='y', alpha=0.5)

    ax.set_ylabel('Production Rate (10$^{24}$ mol s$^{-1}$)', fontsize=11, labelpad=8)

    # Bottom axis
    # Dates x axis
    xlocs = [2458757, 2458788, 2458818, 2458849, 2458880, 2458909, 2458940] # 1st day of month Dec to April
    xlocs = [2458825, 2458825+20, 2458825+40, 2458825+60, 2458825+80, 2458825+100, 2458825+120] # Days post perihelion starting at 0
    xlabs = ['0', '20', '40', '60', '80', '100', '120']
    ax.set(xticks=xlocs, xticklabels=xlabs, xlim=(xlocs[0], xlocs[-1]))
    if j == 2:
        ax.set_xlabel('Days Post Perihelion', fontsize=12, labelpad=10)
    elif j==0 or j==1:
        ax.set(xticklabels=['','','','','','',''], xlim=(xlocs[0], xlocs[-1]))
        ax.tick_params(length=0) # make bottom ticks invisible
    if j == 3:
        xlocs = [2458825-40, 2458825-20, 2458825, 2458825+20, 2458825+40, 2458825+60, 2458825+80, 2458825+100, 2458825+120] # Days post perihelion starting at 0
        xlabs = ['-40', '-20', '0', '20', '40', '60', '80', '100', '120']
        ax.set(xticks=xlocs, xticklabels=xlabs, xlim=(xlocs[0], xlocs[-1]))
        ax.set_ylabel("S' (%/1000$\AA{}$)", fontsize=12, labelpad=8)
        ax.set_xlabel('Days to Perihelion', fontsize=12, labelpad=10)


    # Top axis
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    new_tick_locations = [2458825,2458855,2458867.5,2458886,2458901,2458915,2458928,2458940] # JD corresponding to the Rh values in new_tick_labels
    new_tick_labels = ['0.0', '0.1', '0.2', '0.4', '0.6', '0.8', '1.0', '1.2']
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(new_tick_labels)
    if j == 0:
        ax2.set_xlabel('Distance from Perihelion (au)', fontsize=12, labelpad=10)
    elif j==1 or j==2:
        ax2.set(xticklabels=['','','','','','','',''], xlim=(xlocs[0], xlocs[-1]))
        ax2.tick_params(length=0)
    if j==3:
        new_tick_locations = [2458797,2458825,2458855,2458867.5,2458886,2458901,2458915,2458928] # JD corresponding to the Rh values in new_tick_labels
        new_tick_labels = [ '-0.1', '0.0', '0.1', '0.2','0.4','0.6','0.8','1.0']
        ax2.set_xticks(new_tick_locations)
        ax2.set_xticklabels(new_tick_labels)
        ax2.set_xlabel('Distance to Perihelion (au)', fontsize=12, labelpad=10)
    #ax2.tick_params(length=8)

    

#fig.tight_layout()
bbox = (fig.bbox_inches.from_bounds(-0.1,1,6.5,13)) 
plt.savefig('../figures/ChangeInProdrates.pdf', bbox_inches = bbox)
