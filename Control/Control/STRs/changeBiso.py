# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 19:28:14 2021

@author: Andrew Crossman
"""
#import pandas as pd
#import numpy as np
#import pylab as plt
#import scipy as sp
import re

infile = "Sphere5050_AuNi_STR.stru"
outfile = "Sphere5050_AuNi_STR_Biso.stru"

e1 = 'AU'
e2 = 'NI'
e3 = 'PT'

# Biso values are from

b1 = '0.592500,'
b2 = '0.340100,'
b3 = '0.355700,'

ifile = open(infile, 'r')
ofile = open(outfile, 'w')

data = ifile.readlines()
ofile.write(data[0])
ofile.write(data[1])
ofile.write(data[2])
ofile.write(data[3])

for line in data[4:]:
    l = re.split(r'(\s+)', line)
    if l[0] == e1:
        text = l[0]+l[1]+l[2]+l[3]+l[4]+l[5]+l[6]+l[7]+b1+l[9]+l[10]+l[11]+l[12]+l[13]+l[14]+l[15]+l[16]+l[17]+l[18]+l[19]+l[20]+l[21]+l[22]+l[23]+l[24]+'\n'
        ofile.write(text)
    elif l[0] == e2:
        text = l[0]+l[1]+l[2]+l[3]+l[4]+l[5]+l[6]+l[7]+b2+l[9]+l[10]+l[11]+l[12]+l[13]+l[14]+l[15]+l[16]+l[17]+l[18]+l[19]+l[20]+l[21]+l[22]+l[23]+l[24]+'\n'
        ofile.write(text)
    elif l[0] == e3:
        text = l[0]+l[1]+l[2]+l[3]+l[4]+l[5]+l[6]+l[7]+b3+l[9]+l[10]+l[11]+l[12]+l[13]+l[14]+l[15]+l[16]+l[17]+l[18]+l[19]+l[20]+l[21]+l[22]+l[23]+l[24]+'\n'
        ofile.write(text)

ifile.close()
ofile.close()
