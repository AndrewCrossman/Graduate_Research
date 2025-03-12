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
RMC_runs_path = "C:\\Users\\crossman\\Desktop\\Reverse_Monte_Carlo\\Tests\\Au50B50_CoreShells\\Coreshell5050_PtAu_onAu\\"
RMC_runs = glob.glob(RMC_runs_path+"*1-20-25_v*[!Results]") 
print(RMC_runs_path.split('\\')[-2])
#%%
'''
Extract MMC (Target) Data
'''
MMC_data = Functions.get_stru(RMC_runs_path+"MMC_5000\\2000_voidless.stru",100)
MMC_x, MMC_y = Functions.get_Comp_by_atom_and_abs_dist(MMC_data,num_bins,atoms_per_bin)

#%%
'''
Extract Each RMC (model) Data. The reduced radial distance is used as the x-axis
and the AU Fraction of the Model is used as the y-axis. The x-axis of each model 
is interpolated onto the same scheme as the Target (i.e. the binning is slightly
different run to run so this will put eveything on the same bins).
This effectively stretches the RMC data to larger values of the reduced distance
'''
All_data = []
for path in RMC_runs:
    RMC_data = Functions.get_stru(path+"\\1-20.stru",99)
    # Find Bins by atom count and reduced distance
    x,y = Functions.get_Comp_by_atom_and_abs_dist(RMC_data,num_bins,atoms_per_bin)
    # Interpolates RMC data onto the same x axis as the MMC data
    All_data.append(np.interp(MMC_x, x,y))
# Find the mean/average composition profile from all the RMC runs 
RMC_avg = np.mean(All_data,axis=0)

#%%
'''
Plots the Composition Profiles
'''
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

#%%
'''
Create a csv with the following column headers: reduced_radial_distance, Target_Composition,
Avgerage_RMC_Composition, Run1_Composition, Run2_Composition, etc.
'''
csv = pd.DataFrame({"red_dist": MMC_x, "Target_Comp": MMC_y, "Avg_RMC_Comp": RMC_avg})
csv[["Run{}".format(i+1) for i in range(len(All_data))]] = np.asarray(All_data).T
csv.to_csv(RMC_runs_path+"Composition_Profile_Data_{}.csv".format(RMC_runs_path.split('\\')[-2]),
           sep='\t', index=False)
#%%
'''
Finds the sum of the squared difference between the Average_RMC_Compostion and the Target_Composition
'''
difference = np.subtract(np.array(MMC_y),np.array(RMC_avg))
difference_squared = np.square(difference)
print('X^2: {}'.format(sum(difference_squared)))

R_wp = np.sqrt(sum(difference*difference)/sum(np.array(MMC_y)*np.array(MMC_y)))*100
print("R_wp: {}%".format(round(R_wp,4)))