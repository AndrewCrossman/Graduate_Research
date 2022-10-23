#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 20:13:39 2022

@author: andrew
"""
import matplotlib.pyplot as plt
import numpy as np

#ext = "Tests/Sigma_0070_2times200_Homo_exact/"
ext=""
infile = "suite.log"
#infile = ext + "suite.log"

ifile = open(infile, 'r')

lines = ifile.readlines()#[85:]
data = []

for line in lines:
    l = line.split()
    print(l)
    if l==[]:
        pass
    elif l[0]=="Lennard":
        data.append(float(l[3])*float(l[4])+float(l[6])*float(l[7]))
        '''if l[3]=="********":
            
        else:
            data.append(float(l[3])*float(l[4])+float(l[6])*float(l[7]))'''
        
print(min(data))
fig, ax = plt.subplots(figsize=(9,6))
ax.scatter(range(len(data)), data, color='blue', marker='o',label=r'E') 
ax.set(xlabel='MC Cycles', ylabel=r'$\delta E$',
       title='Change in '+r'$\delta E $'+' MC vs Cycles')
ax.legend(loc="upper right")
ax.grid()
fig.savefig(ext+"energy.pdf",dpi=100)
plt.tight_layout()
plt.show()
#plt.save("energy.png")
