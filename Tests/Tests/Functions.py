# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 17:11:16 2024
@author: crossman
"""
#%%
'''
IMPORTS
'''
import shutil                   #allows uploading files to drive
import glob                     #useful for collecting and organizing files
import csv
from pathlib import Path        #useful for reading text from files
import numpy as np              #math library
import pandas as pd             #data organization library
import scipy as sp              #math library
from scipy import interpolate
import re                       #string processing library
from scipy import optimize      #math fitting library
import matplotlib.pylab as plt  #allows for plotting
#%%
'''
Reads in DISCUS stru file and formats it into a pandas dataframe. Data is then 
processed. Percentile is used as the quantity to which to set the reduced distance to 1
In other words, 98 would mean that a reduced distance of 1 will be at the 98th
percentile. This helps to account for trace amounts of evaporating atoms from the 
RMC reconstructions. Returns processed dataframe.
'''
def get_stru(filename,percentile):
    data = pd.read_csv(filename, skiprows=4, sep="\\s+")
    data = data.drop(data.columns[4:], axis=1)
    data = data.replace(',','', regex=True)
    data = data.set_axis(["atom", "x", "y", "z"], axis=1)
    data = data[data.atom!='VOID']
    data["x"] = data["x"].astype(float)
    data["y"] = data["y"].astype(float)
    data["z"] = data["z"].astype(float)
    data["distance"] = np.sqrt((data["x"]*data["x"] + data["y"]*data["y"] + data["z"]*data["z"])*(4.0789*4.0789))
    data["reduced_distance"] = data["distance"]*(1/np.percentile(data["distance"], percentile))
    # Output relevant structural info
    print("Min-Max distance in angstroms: {} - {}".format(str(round(data["distance"].min(),5)),str(round(data["distance"].max(),5))))
    return data

#%%
'''
Reads in a processed pandas dataframe and finds the composition profile. The number
of bins to use is numb_bins and the number of atoms per bin is atoms_per_bin.
Returns the composition (y) at distance (x)
'''
def get_Comp_by_atom_and_abs_dist(dataframe, num_bins, atoms_per_bin):
    x, y = [], []
    for i in range(num_bins):
        # Finds the desired atoms 'atoms_per_bin'
        RMC_temp_red = dataframe.nsmallest(atoms_per_bin*(i+1), "reduced_distance")[-atoms_per_bin:]
        # Solves for the mean reduced distance for all atoms belonging to the the bin
        x.append(RMC_temp_red["reduced_distance"].mean())
        # Solves for the Au fraction of the bin
        y.append(len(RMC_temp_red[RMC_temp_red["atom"]=='AU']) / len(RMC_temp_red))
    return x,y
#%%
'''
Reads in a file of text and finds every chi_squared value. Outputs an array of 
the chi_squared values in chronological order.
'''
def get_Chi_Squared(filename):
    file_text = Path(filename).read_text()
    str_to_search = "s2x2"
    all_X2s = np.array(re.findall("s2x2:\ *(\d+\.?\d*[Ee]?-?\d+)", file_text),dtype=float)
    return all_X2s
#%%

