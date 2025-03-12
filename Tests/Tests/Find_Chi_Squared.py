# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 10:01:06 2024
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
from astropy.coordinates import cartesian_to_spherical
import re                       #string processing library
from scipy import optimize      #math fitting library
import matplotlib.pylab as plt  #allows for plotting
# Script Imports
'''
In the upper right on anaconda move working directory to where this file and Functions.py are located
'''
import Functions
#%%
'''
Important Variables Go Here
'''
#%%
'''
Find RMC Run Directories
'''
RMC_runs_path = "C:\\Users\\crossman\\Desktop\Reverse_Monte_Carlo\\Tests\\Au50B50_CoreShells\Coreshell5050_AuNi_onAu\\"
RMC_runs = glob.glob(RMC_runs_path+"*1-20-25_v*[!Results]")    # [!Results] excludes Results in the string when searching
print(RMC_runs)
#%%
'''
Extract all relevant chi squared data from log file
'''
fig1, ax1 = plt.subplots(figsize=(6,4,),dpi=200)
left, bottom, width, height = [0.5, 0.5, 0.45, 0.35]
ax2 = fig1.add_axes([left, bottom, width, height])
Final_X2s = []
for path in RMC_runs:
    RMC_data = Functions.get_Chi_Squared(path+"\\rmc_suite.log")
    print(len(RMC_data))
    ax1.plot([x *25/4 for x in range(len(RMC_data))],RMC_data, color='blue', alpha=0.5)
    ax2.plot([x *25/4 for x in range(len(RMC_data))][20:60],RMC_data[20:60], color='blue', alpha=0.5)
    Final_X2s.append(RMC_data[-1])
print(np.mean(Final_X2s))
# Plot settings
ax1.set_title("{} ".format(RMC_runs_path.split('\\')[-2])+r'$\chi^2$ Evolution', fontsize=20)
ax1.set_ylabel(r'$\chi^2$', fontsize=16)
ax1.tick_params(axis='y', labelsize=12, direction='in', length=8)
ax1.set_yscale('log')
ax1.set_xlabel(r'RMC Cycles', fontsize=16)
ax1.tick_params(axis='x', labelsize=12, direction='in', length=8)
ax1.grid(True)
ax1.set_axisbelow(True)
ax2.grid(True)
ax2.set_axisbelow(True)
plt.tight_layout()