#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 20:13:39 2022

@author: andrew
"""
import matplotlib.pyplot as plt
import numpy as np

#ext = "Tests/Sigma_0070_2times200_Homo_exact/"
ext="MC_Coreshell5050_PtAu_onAu_r25_5000Cycles/"
infile = ext+"suite.log"
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
print(len(data))
x = [1000*(x+1) for x in range(len(data))]
fig, ax = plt.subplots(figsize=(9,6))
ax.scatter(x, data, color='blue', marker='o',label=r'E') 
ax.set(xlabel='MC Cycles', ylabel=r'Energy (a.u.)',
       title='PtAu: Total Lattice Energy'+' vs MC Cycles')
ax.legend(loc="upper right")
ax.grid()
fig.savefig("PtAu_r25_energy.pdf",dpi=100)
plt.tight_layout()
plt.show()
#plt.save("energy.png")
