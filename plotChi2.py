#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  1 19:56:14 2022

@author: andrew
"""
import matplotlib.pyplot as plt
import numpy as np

ext = "RMC_Coreshell5050_PtAu_on_Au_1-10-10/"
infile = ext + "rmc_suite.log"

ifile = open(infile, 'r')

lines = ifile.readlines()[85:]
data = []

for line in lines:
    l = line.split()
   # print(l)
    if l==[]:
        pass
    elif l[0]=="Gen:":
        data.append(float(l[10]))
        
print(min(data))
print(len(data))
fig, ax = plt.subplots(figsize=(9,6))
ax.scatter(range(len(data)), data, color='blue', marker='o',label=r'$\chi ^2$') 
ax.set(xlabel='Cycles', ylabel=r'$\chi ^2$',
       title='PtAu: Change in '+r'$\chi ^2$'+' vs Cycles')
ax.legend(loc="upper right")
ax.grid()
fig.savefig(ext+"PtAu_chi2.pdf",dpi=100)
plt.tight_layout()
plt.show()
        

