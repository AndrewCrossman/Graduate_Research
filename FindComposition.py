#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 12:40:49 2022

@author: andrew
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

numbin = 15     #number of bins/slices of sphere
lp = 4.0789
ext =""
#file="Tests/RMC_Coreshell5050_AuPd_on_Au_1-10-10/1-10.stru"
file="RMC_Coreshell5050_PtAu_on_Au_1-10-10/1-1.stru"
title = "Composition of RMC PtAu (10Cycles)"
savename = "PtAu_CompProf_1-1.pdf"


data = pd.read_csv(file, skiprows=3)
data = data.drop(data.columns[3:], axis=1)
data[["atoms", "x"]] = data[data.columns[0]].str.split(expand=True)
data = data.drop(data.columns[0], axis=1)
data = data.set_axis(["y", "z", "atom", "x"], axis=1, inplace=False)
data["x"] = data["x"].astype(float)
data["y"] = data["y"].astype(float)
data["z"] = data["z"].astype(float)
data["distance"] = np.sqrt((data["x"]*data["x"] + data["y"]*data["y"] + data["z"]*data["z"])*(lp*lp))
print("Max distance in angstroms "+str(data["distance"].max()))
print("Min distance in angstroms "+str(data["distance"].min()))

realBins = np.linspace(0,data["distance"].max()+.01,numbin)
data["bin"] = pd.cut(data["distance"], bins=realBins, right=False, precision=2,include_lowest='True')

data2 = data.groupby("bin",observed=True)["atom"].value_counts(normalize=True).unstack()

tracker = 0
colors = ['red','blue']
fig, ax = plt.subplots(figsize=(6,4))
for key in data2.keys():
    data3 = data2[key].reset_index()
    ax.scatter(range(len(data3['bin'])),data3[key],color=colors[tracker])
    tracker+=1
x = []
for c in data3.bin:
    x.append(c)

realBins_truncated = np.around(realBins,decimals=2)
ax.set_title(title, fontsize=20)
ax.set_ylabel("Fraction", fontsize=16)
ax.set_xlim(-.5,30.5)
ax.set_ylim(-.05,1.05)
ax.tick_params(axis='y', labelsize=12)
ax.set_xticks(np.arange(len(x)))
ax.set_xticklabels(x)
ax.set_xticklabels(realBins_truncated[1:])
ax.set_xlabel("Distance [A]", fontsize=16)
ax.tick_params(axis='x', labelsize=10, labelrotation = 30)
ax.legend()
ax.grid(True)
plt.tight_layout()
plt.savefig(savename, dpi=400)



