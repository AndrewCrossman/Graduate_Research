# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 12:26:16 2024
@author: crossman
"""
#%%
'''
IMPORTS
'''
# Library Imports
import os                       #useful for importing other scripts
import shutil                   #allows uploading files to drive
import glob                     #useful for collecting and organizing files
import csv
import numpy as np              #math library
import pandas as pd             #data organization library
import scipy as sp              #math library
from scipy import interpolate
import re                       #string processing library
from scipy import optimize      #math fitting library
import matplotlib.pylab as plt  #allows for plotting
# Script Imports
import Functions
#%%
'''
Important Variables Go Here
'''
num_bins, atoms_per_bin = 25,157
#%%
'''
Find RMC Run Directories
'''
RMC_runs_path = "C:\\Users\\crossman\\Desktop\Reverse_Monte_Carlo\\Tests\\Au50Pd50_Profiles\\Dt500\\"
RMC_runs = glob.glob(RMC_runs_path+"*1-20-25_v*[!Results]") 
print(RMC_runs_path.split('\\')[-2])
#%%
'''
Extract MMC (Target) Data
'''
MMC_data = Functions.get_df(RMC_runs_path+"MMC_2000\\2000_voidless.stru",100)
MMC_x, MMC_y = Functions.get_Comp_by_atom_and_abs_dist(MMC_data,num_bins,atoms_per_bin)

#%%
All_data = []
for path in RMC_runs:
    RMC_data = Functions.get_df(path+"\\1-20.stru",98)
    # Find Bins by atom count and reduced distance
    x,y = Functions.get_Comp_by_atom_and_abs_dist(RMC_data,num_bins,atoms_per_bin)
    All_data.append(np.interp(MMC_x, x,y))
RMC_avg = np.mean(All_data,axis=0)
    
#%%

fig1, ax1 = plt.subplots(figsize=(6,4,),dpi=200)
for y in All_data:
    ax1.plot(MMC_x,y, color="blue",alpha=.5, zorder=0)
ax1.plot(MMC_x,MMC_y, color='black', label = "Target", zorder=1)
ax1.plot(MMC_x,RMC_avg, color='red', label = "Model", zorder=1)
# Plot settings
ax1.set_title("{} By Atom Count".format(RMC_runs_path.split('\\')[-2]), fontsize=20)
ax1.set_ylabel("Au Fraction", fontsize=16)
ax1.tick_params(axis='y', labelsize=12, direction='in', length=8)
ax1.set_xlabel(r'Reduced Distance ($r/r_0$)', fontsize=16)
ax1.tick_params(axis='x', labelsize=12, direction='in', length=8)
ax1.legend()
ax1.grid(True)
ax1.set_axisbelow(True)
plt.tight_layout()