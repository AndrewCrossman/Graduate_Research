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
RMC_runs_path = "C:\\Users\\crossman\\Desktop\Reverse_Monte_Carlo\\Tests\\Au50Pd50_Profiles\\Dt2\\"
RMC_runs = glob.glob(RMC_runs_path+"*1-20-25_v*[!Results]") 
print(RMC_runs_path.split('\\')[-2])
#%%
'''
Extract MMC (Target) Data
'''
MMC_pdf = pd.read_csv(RMC_runs_path + "MMC_2000\\2000.apd", names=['r','G','idk1','idk2'], header=None, sep="\\s+")
MMC_xrd = pd.read_csv(RMC_runs_path + "MMC_2000\\2000.xrd", names=['tth','I'], header=None, sep="\\s+")
MMC_bld11 = pd.read_csv(RMC_runs_path + "MMC_2000\\2000_AuAu.bld", names=['r','y'], header=None, sep="\\s+") 
MMC_bld12 = pd.read_csv(RMC_runs_path + "MMC_2000\\2000_AuPd.bld", names=['r','y'], header=None, sep="\\s+") 
MMC_bld22 = pd.read_csv(RMC_runs_path + "MMC_2000\\2000_PdPd.bld", names=['r','y'], header=None, sep="\\s+") 
#%%
'''
Find Average PDF, XRD, BLDs of RMC runs
'''
RMC_pdf = pd.DataFrame()
RMC_xrd = pd.DataFrame()
RMC_bld11 = pd.DataFrame()
RMC_bld12 = pd.DataFrame()
RMC_bld22 = pd.DataFrame()
for path in RMC_runs:
    pdf = pd.read_csv(RMC_runs_path + "1-20.apd", names=['r','G','idk1','idk2'], header=None, sep="\\s+")
    xrd = pd.read_csv(RMC_runs_path + "1-20.xrd", names=['tth','I'], header=None, sep="\\s+")
    bld11 = pd.read_csv(RMC_runs_path + "1-20_AuAu.bld", names=['r','y'], header=None, sep="\\s+")
    bld12 = pd.read_csv(RMC_runs_path + "1-20_AuPd.bld", names=['r','y'], header=None, sep="\\s+")
    bld22 = pd.read_csv(RMC_runs_path + "1-20_PdPd.bld", names=['r','y'], header=None, sep="\\s+")
    
