# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 14:35:18 2024
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
RMC_runs_path = "C:\\Users\\crossman\\Desktop\\Reverse_Monte_Carlo\\Tests\\Au50Pd50_Profiles\\Dt2\\"
RMC_runs = glob.glob(RMC_runs_path+"*1-20-25_v*[!Results]") 
print(RMC_runs_path.split('\\')[-2])
#%%
'''
Extract MMC (Target) Data
'''
fig, axs = plt.subplots(6,dpi=400,figsize=(7,10))
i=0
MMC_data,MMC_r0 = Functions.get_stru(RMC_runs_path+"MMC_2000\\2000_voidless.stru",100,True)
MMC_x, MMC_y = Functions.get_Comp_by_atom_and_abs_dist(MMC_data,num_bins,atoms_per_bin)
axs[i].hist(MMC_data["distance"],bins=40, range=(-1, 30))

i=i+1
#%%
'''
Plots Binned Data
'''
All_data = []
r_0s = []
for path in RMC_runs:
    RMC_data,r_0 = Functions.get_stru(path+"\\1-20.stru",100,True)
    axs[i].hist(RMC_data["distance"],bins=40, range=(-1, 30))
    i=i+1
#%%
#%%
'''
Extract MMC (Target) Data
'''
MMC_data,MMC_r0 = Functions.get_stru(RMC_runs_path+"MMC_2000\\2000_voidless.stru",100,True)
MMC_x, MMC_y = Functions.get_Comp_by_atom_and_abs_dist(MMC_data,num_bins,atoms_per_bin)



#%%
'''
Extract Each RMC (model) Data. The reduced radial distance is used as the x-axis
and the AU Fraction of the Model is used as the y-axis. The x-axis of each model 
is interpolated onto the same scheme as the Target (i.e. the binning is slightly
different run to run so this will put eveything on the same bins).
This effectively stretches the RMC data to larger values of the reduced distance
'''
print('############################')
print('Analyzing Composition Profile')
All_data = []
r_0s = []
for path in RMC_runs:
    RMC_data,r_0 = Functions.get_stru(path+"\\1-20.stru",100,True)
    # Find Bins by atom count and reduced distance
    x,y = Functions.get_Comp_by_atom_and_abs_dist(RMC_data,num_bins,atoms_per_bin)
    #x,y = Functions.get_Comp_by_atom_and_abs_dist(RMC_data,num_bins,atoms_per_bin,atom='CU')
    # Interpolates RMC data onto the same x axis as the MMC data
    All_data.append(np.interp(MMC_x, x,y))
    r_0s.append(r_0)
# Find the mean/average composition profile from all the RMC runs 
print('Checking 100 percentile...')
best_RMC_avg = np.mean(All_data,axis=0)
best_R_wp = Functions.get_R_wp(MMC_y, best_RMC_avg)
best_percent=100
best_All_data = All_data
best_r_0_avg, r_0_original = np.mean(r_0s), np.mean(r_0s)
#%%
'''
Modify the reduced radial distance to find the best effective radial distance 
that fits the composition plot. Essentially this does the exact same thing as the
cell above but identifies the best percentile to use.
'''
for i in np.arange(90,100,.5):
    print('Checking {} percentile...'.format(i))
    All_data = []
    r_0s = []
    for path in RMC_runs:
        RMC_data, r_0 = Functions.get_stru(path+"\\1-20.stru",i,False)
        # Find Bins by atom count and reduced distance
        x,y = Functions.get_Comp_by_atom_and_abs_dist(RMC_data,num_bins,atoms_per_bin)
        # Interpolates RMC data onto the same x axis as the MMC data
        All_data.append(np.interp(MMC_x, x,y))
        r_0s.append(r_0)
    # Find the mean/average composition profile from all the RMC runs 
    RMC_avg = np.mean(All_data,axis=0)
    r_0_avg = np.mean(r_0s)
    R_wp = Functions.get_R_wp(MMC_y, RMC_avg)
    if R_wp < best_R_wp:
        best_R_wp = R_wp
        best_RMC_avg = RMC_avg
        best_percent = i
        best_r_0_avg = r_0_avg
        best_All_data = All_data
print('############################')
print('MMC_r_0={} - Best r_0={}/{}=.{}, R_wp={}%'.format(round(MMC_r0,4),
      round(best_r_0_avg,4),round(r_0_original,4),best_percent,round(best_R_wp,4)))
print('############################')