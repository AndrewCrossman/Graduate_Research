# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 23:53:21 2024
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
Find RMC Run Directories
'''
RMC_runs_path = "C:\\Users\\crossman\\Desktop\\Reverse_Monte_Carlo\\Tests\\Au50B50_CoreShells\\Coreshell5050_PtAu_onAu\\"
RMC_runs = glob.glob(RMC_runs_path+"*1-20-25_v*[!Results]") 
#print(RMC_runs)
#print(RMC_runs_path.split('\\')[-2])
#%%
'''
Extract MMC (Target) Data
'''
MMC_pdf = pd.read_csv(RMC_runs_path + "MMC_5000\\2000.apd", names=['r','G','idk1','idk2'], header=None, sep="\\s+")
MMC_xrd = pd.read_csv(RMC_runs_path + "MMC_5000\\2000.xrd", names=['tth','I'], header=None, sep="\\s+")
MMC_bld11 = pd.read_csv(RMC_runs_path + "MMC_5000\\2000_AuAu.bld", names=['r','y'], header=None, sep="\\s+") 
MMC_bld12 = pd.read_csv(RMC_runs_path + "MMC_5000\\2000_AuPt.bld", names=['r','y'], header=None, sep="\\s+") 
MMC_bld22 = pd.read_csv(RMC_runs_path + "MMC_5000\\2000_PtPt.bld", names=['r','y'], header=None, sep="\\s+") 
#%%
'''
Find Average PDF, XRD, BLDs of RMC runs
'''
# create five arrays to hold measurements from the different runs
RMC_pdf = []
RMC_xrd = []
RMC_bld11 = []
RMC_bld12 = []
RMC_bld22 = []
for path in RMC_runs:
    # extract data for the current 'path'
    pdf_data = pd.read_csv(path + "\\1-20.apd", names=['r','G','idk1','idk2'], header=None, sep="\\s+")
    xrd_data = pd.read_csv(path + "\\1-20.xrd", names=['tth','I'], header=None, sep="\\s+")
    bld11_data = pd.read_csv(path + "\\1-20_AuAu.bld", names=['r','y'], header=None, sep="\\s+")
    bld12_data = pd.read_csv(path + "\\1-20_AuPt.bld", names=['r','y'], header=None, sep="\\s+")
    bld22_data = pd.read_csv(path + "\\1-20_PtPt.bld", names=['r','y'], header=None, sep="\\s+")
    # append extracted data to arrays that hold all of the run-to-run data after interpolating onto 
    # the x-axis used for the MMC. This way all values are on the same x-axis spacings
    RMC_pdf.append( np.interp(MMC_pdf['r'],pdf_data['r'],pdf_data['G']) )
    RMC_xrd.append( np.interp(MMC_xrd['tth'],xrd_data['tth'],xrd_data['I']) )
    RMC_bld11.append( np.interp(MMC_bld11['r'],bld11_data['r'],bld11_data['y']) )
    RMC_bld12.append( np.interp(MMC_bld12['r'],bld12_data['r'],bld12_data['y']) )
    RMC_bld22.append( np.interp(MMC_bld22['r'],bld22_data['r'],bld22_data['y']) )
# calculate the average PDF, XRD, BLDs by averaging on a point-by-point basis
RMC_pdf_avg = np.mean(RMC_pdf, axis=0)
RMC_xrd_avg = np.mean(RMC_xrd, axis=0)
RMC_bld11_avg = np.mean(RMC_bld11, axis=0)
RMC_bld12_avg = np.mean(RMC_bld12, axis=0)
RMC_bld22_avg = np.mean(RMC_bld22, axis=0)
#%%
'''
Find R_wp values (good fit if <5%)
'''
pdf_dif = MMC_pdf['G'] - RMC_pdf_avg
xrd_dif = MMC_xrd['I'] - RMC_xrd_avg
bld11_dif = MMC_bld11['y'] - RMC_bld11_avg
bld12_dif = MMC_bld12['y'] - RMC_bld12_avg
bld22_dif = MMC_bld22['y'] - RMC_bld22_avg
# Calculate R_wp values
R_wp_pdf = np.sqrt(np.sum(pdf_dif*pdf_dif)/np.sum(MMC_pdf['G'] * MMC_pdf['G'])) *100
R_wp_xrd = np.sqrt(np.sum(xrd_dif*xrd_dif)/np.sum(MMC_xrd['I'] * MMC_xrd['I'])) *100
R_wp_bld11 = np.sqrt(np.sum(bld11_dif*bld11_dif)/np.sum(MMC_bld11['y'] * MMC_bld11['y'])) *100
R_wp_bld12 = np.sqrt(np.sum(bld12_dif*bld12_dif)/np.sum(MMC_bld12['y'] * MMC_bld12['y'])) *100
R_wp_bld22 = np.sqrt(np.sum(bld22_dif*bld22_dif)/np.sum(MMC_bld22['y'] * MMC_bld22['y'])) *100
# print R_wp values
print(RMC_runs_path.split('\\')[-2])
print("PDF R_wp: {}%".format(round(R_wp_pdf,4)))
print("XRD R_wp: {}%".format(round(R_wp_xrd,4)))
print("BLD11 R_wp: {}%".format(round(R_wp_bld11,4)))
print("BLD12 R_wp: {}%".format(round(R_wp_bld12,4)))
print("BLD22 R_wp: {}%".format(round(R_wp_bld22,4)))
print('{} {} {} {} {}'.format( round(R_wp_pdf,4), round(R_wp_xrd,4), round(R_wp_bld11,4), round(R_wp_bld12,4), round(R_wp_bld22,4)))
#%%
'''
Plot RMC data vs PDF Data
'''
fig, axs = plt.subplots(5,dpi=400,figsize=(6,8))

axs[0].plot(MMC_pdf['r'],MMC_pdf['G'])
axs[0].plot(MMC_pdf['r'],RMC_pdf_avg)

axs[1].plot(MMC_xrd['tth'],MMC_xrd['I'])
axs[1].plot(MMC_xrd['tth'],RMC_xrd_avg)

axs[2].plot(MMC_bld11['r'],MMC_bld11['y'])
axs[2].plot(MMC_bld11['r'],RMC_bld11_avg)

axs[3].plot(MMC_bld12['r'],MMC_bld12['y'])
axs[3].plot(MMC_bld12['r'],RMC_bld12_avg)
    
axs[4].plot(MMC_bld22['r'],MMC_bld22['y'])
axs[4].plot(MMC_bld22['r'],RMC_bld22_avg)
#%%
