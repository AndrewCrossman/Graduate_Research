# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 16:24:24 2025
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
from sklearn.decomposition import PCA              #used for PCA
from sklearn.preprocessing import StandardScaler   #used for PCA
import matplotlib.pylab as plt  #allows for plotting
from matplotlib.pyplot import cm
# Script Imports
import Functions
#%%
'''
Find directory and list all pdfs of interest
'''
# directory path 
directory_path = "C:\\Users\\crossman\\Desktop\\Reverse_Monte_Carlo\\Tests\Au50Pd50_Profiles\\Dt0\\"
# pdfs of interest
pdfs = glob.glob(directory_path+"RMC_1-20-25_v1\\*.apd")
# ansatz pdf
ansatz_pdf = directory_path+"Sphere5050_AuPd_onAu_r25_voidless.apd"
#%%
'''
Extract and Format Data. Data must assume an [m x n] matrix, where the rows are
independent samples and the columns are features corresponding to each sample
'''
# initialize dataframe
all_pdfs = pd.read_csv(ansatz_pdf, sep="\\s+", usecols=[0,1],
                       header=None, names=['r',0])
r_values = all_pdfs['r']
all_pdfs = all_pdfs.drop('r', axis=1)
# populate dataframe with pdf data
for file in pdfs:
    pdf = pd.read_csv(file, sep="\\s+", usecols=[1], header=None)
    pdf.columns = [int(re.findall("r*\d+-(\d+).apd", file)[0])]
    all_pdfs = pd.concat([all_pdfs, pdf], axis=1)
    
# sort columns by value except for the first column which lists the interatomic distance
all_pdfs = pd.concat([all_pdfs[[all_pdfs.columns[0]]],all_pdfs.reindex(sorted(all_pdfs.columns[1:]), axis=1)],axis=1)
# convert header names to strings
all_pdfs.columns = all_pdfs.columns.astype(str) 
# transpose data so that row headers are the different samples (diffusive states)
# and the columns are the different measurements (interatomic distances)
all_pdfs = all_pdfs.transpose()
#%%
'''
Execute PCA Measurement and print eigenvalue info plus relative explained variance

transformed_data is a 2D array of [sample x PC] where each element is the weight of 
the PC for that sample (i think!)
'''
# Get scaler object, pca object, scaled data (normalized data), and transformed data (PC matrix)
scaler, pca, scaled_data, transformed_data = Functions.get_PCA(all_pdfs,10)
# transformed data in pandas dataframe objext
df_transformed_data = pd.DataFrame(transformed_data)
# The amount of variance explained by each of the selected components (i.e. eigenvalue)
print("Component Eigenvalues:")
print(pca.explained_variance_)
# Percentage of variance explained by each of the selected components in order i.e. PC1 %, PC2 %, ...
print("Component relative explained variance as a % of total:")
print(pca.explained_variance_ratio_)
#eigenvectors = pd.DataFrame(pca.components_).transpose()
#%%
'''
Create a scree plot. This shows the relative explained variance as a % of total (100%).
The x-axis describes the principle component number and the y-axis its value.
'''
# Gets the principle numbers
PC_values = np.arange(pca.n_components_) + 1
fig1, ax1 = plt.subplots(figsize=(6,4,),dpi=400)
ax1.plot(PC_values, pca.explained_variance_ratio_, 'o-', linewidth=2, color='black')
ax1.set_title('Scree Plot', fontsize=14)
ax1.set_xlabel('Principal Component', fontsize=14)
ax1.set_ylabel('Variance Explained', fontsize=14)
ax1.tick_params(axis='y', labelsize=12, direction='in', length=8)
ax1.tick_params(axis='x', labelsize=12, direction='in', length=8)
ax1.grid(True)
ax1.set_axisbelow(True)
plt.tight_layout()
#%%
'''
Creats a plot that shows the cumulative explained variance i.e. what percentage 
of all data can be explained using all components up to the amount along the x-axis.
For example, if at x=3, y=.97 then using components 1,2,3 97% of data can be calulcated
'''
fig2, ax2 = plt.subplots(figsize=(6,4,),dpi=400)
ax2.plot(PC_values, np.cumsum(pca.explained_variance_ratio_), 'o-', linewidth=2, color='black')
ax2.set_title('Cumulative Explained Variance')
ax2.set_xlabel('Principle Components', fontsize=14)
ax2.set_ylabel('Cumulative Explained Variance', fontsize=14);
ax2.tick_params(axis='y', labelsize=12, direction='in', length=8)
ax2.tick_params(axis='x', labelsize=12, direction='in', length=8)
ax2.grid(True)
ax2.set_axisbelow(True)
plt.tight_layout()
#%%
'''
Creates a scatter plot of all the principle components against one another. This
is useful for identifying correlations between different components. If uncorrelated,
the plot should look unstructure or uniform. Otherwise, there is likely some sort
of correlation. The diagonal elements should look like gaussians if uncorrelated. 
'''
#pd.plotting.scatter_matrix(pd.DataFrame(transformed_data), color='black', alpha=.5, hist_kwds={'color':'black'});
fig3, axes3 = plt.subplots(PC_values[-1], PC_values[-1], figsize=(12,12), dpi=400)
for i in range(len(axes3)):
    for j in range(len(axes3[0])):
        if i==j:
            axes3[i][j].hist(df_transformed_data[i])            
        else:
            axes3[i][j].scatter(df_transformed_data[i],df_transformed_data[j],
                               c=cm.winter(np.linspace(0, 1, len(transformed_data[:, 1]))),
                               edgecolor='none', alpha=0.5,
                               cmap=plt.cm.get_cmap('winter', len(transformed_data[:, 1])));
        if i==0:
            axes3[i][j].set_title('PC{}'.format(j+1))
        if j==0:
            axes3[i][j].set_ylabel('PC{}'.format(i+1))
plt.tight_layout()
#%%
'''
Creates a plot of the weighting (eigenvector elements) as a function of the
associated samples i.e. the exact weights needed for each eignvector such that:
    scaled_data = a_1*v_1 + a_2*v_2 + ...
where a_i are the weights and v_i are the eigenvectors. (i think!)
'''
fig4, ax4 = plt.subplots(figsize=(6,4,),dpi=400)
cmap = plt.cm.get_cmap('viridis') # You can choose any colormap
num_colors = len(transformed_data[0])  # Number of color divisions for the axes
for i in range(num_colors):
    ax4.plot(transformed_data[:, i], 'o-', color = cmap(i / num_colors), label="PC{}".format(i+1))
ax4.set_title("PC Weights vs Cycles", fontsize=14)
ax4.set_xlabel('Cycles/25', fontsize=14)
ax4.set_ylabel('Weighting', fontsize=14)
ax2.tick_params(axis='y', labelsize=12, direction='in', length=8)
ax2.tick_params(axis='x', labelsize=12, direction='in', length=8)
ax4.legend()
ax4.grid(True)
plt.tight_layout()
plt.show()
#%%
'''
Plot that shows the individual eigenvectors in normalized coordinates and in 
original coordinates. 
'''
eigenvectors = pca.components_
#converts eignevectors from normalized to original units
eigenvectors_og_coords = pd.DataFrame(scaler.inverse_transform(eigenvectors).transpose())

fig9, axes9 = plt.subplots(PC_values[-1], 1, figsize=(12,12), sharex=True, dpi=400);
for i in range(len(axes9)):
    axes9[i].plot(r_values,eigenvectors[i], color = cmap(i / num_colors), label="PC{}".format(i+1));
    axes9[i].set_ylabel(r'$G_{norm}$ (r)', fontsize=14);
    axes9[i].tick_params(axis='y', labelsize=12, direction='in', length=8);
    axes9[i].tick_params(axis='x', labelsize=12, direction='in', length=8);
    axes9[i].grid(True);
    axes9[i].set_axisbelow(True);
    axes9[i].legend(loc='lower right');
    if i==0:
        axes9[i].set_title('Eigenvectors in Normalized Coordinates', fontsize=14);
    if i==len(axes9)-1:
        axes9[i].set_xlabel(r'r ($\AA$)', fontsize=14);
plt.tight_layout();
plt.show();

fig5, axes5 = plt.subplots(PC_values[-1], 1, figsize=(12,12), sharex=True, dpi=400);
for i in range(len(axes5)):
    axes5[i].plot(r_values,eigenvectors_og_coords[i], color = cmap(i / num_colors), label="PC{}".format(i+1));
    axes5[i].set_ylabel(r'G (r)', fontsize=14);
    axes5[i].tick_params(axis='y', labelsize=12, direction='in', length=8);
    axes5[i].tick_params(axis='x', labelsize=12, direction='in', length=8);
    axes5[i].grid(True);
    axes5[i].set_axisbelow(True);
    axes5[i].legend(loc='lower right');
    if i==0:
        axes5[i].set_title('Eigenvectors in Original Coordinates', fontsize=14);
    if i==len(axes5)-1:
        axes5[i].set_xlabel(r'r ($\AA$)', fontsize=14);
plt.tight_layout();
plt.show();

fig6, ax6 = plt.subplots(figsize=(6,4,),dpi=400)
for i in range(len(axes5)):
    ax6.plot(r_values[3540:3550],eigenvectors_og_coords[i][3540:3550], label="PC{}".format(i+1))
ax6.set_title('Eigenvectors in Original Coordinates')
ax6.set_xlabel(r'r ($\AA$)', fontsize=14)
ax6.set_ylabel('G(r)', fontsize=14)
ax6.tick_params(axis='y', labelsize=12, direction='in', length=8)
ax6.tick_params(axis='x', labelsize=12, direction='in', length=8)
ax6.legend()
ax6.grid(True)
ax6.set_axisbelow(True)
plt.tight_layout()
plt.show()
#%%
'''
Plot that compare original data to modeled data based on principle components 
'''
# Gets original unscaled data (pdfs_all)
y1 = all_pdfs.transpose()
# Gets unscaled model data obtained through principle components
model_pdfs = pd.DataFrame(    # convert to dataframe
        scaler.inverse_transform(   # convert to original coordinates
            pd.DataFrame(   # convert to dataframe
                pca.inverse_transform(transformed_data).transpose() # convert principle components to scaled coordinates
            ).transpose() 
        )
     ).transpose()

fig7, ax7 = plt.subplots(figsize=(6,4,),dpi=400)
ax7.plot(r_values, y1.iloc[:, 1], color='black', label='Exp')
ax7.plot(r_values, model_pdfs.iloc[:, 1], color='red', label='PCA')
ax7.set_title('Experimental vs Model G(r)')
ax7.set_xlabel(r'r ($\AA$)', fontsize=14)
ax7.set_ylabel('G(r)', fontsize=14)
ax7.tick_params(axis='y', labelsize=12, direction='in', length=8)
ax7.tick_params(axis='x', labelsize=12, direction='in', length=8)
ax7.legend()
ax7.grid(True)
ax7.set_axisbelow(True)
plt.tight_layout()
plt.show()
#%%
'''
Plot that compares the weighted R-factors (R_wp)
'''
R_wps = []
for i in range(len(model_pdfs.columns)):
    R_wps.append(Functions.get_R_wp(y1.iloc[:, i], model_pdfs.iloc[:, i]))

fig8, ax8 = plt.subplots(figsize=(6,4,),dpi=400)
ax8.plot(R_wps, color='black')
ax8.set_title(r'Experimental vs Model $R_{wp}$')
ax8.set_xlabel(r'Cycle/25', fontsize=14)
ax8.set_ylabel(r'$R_{wp}$ (%)', fontsize=14)
ax8.tick_params(axis='y', labelsize=12, direction='in', length=8)
ax8.tick_params(axis='x', labelsize=12, direction='in', length=8)
ax8.grid(True)
ax8.set_axisbelow(True)
plt.tight_layout()
plt.show()