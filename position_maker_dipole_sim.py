#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 16:06:37 2022

@author: maedeh
"""

import numpy as np

def pos(a): 
    
    y = [-300]
    x = [-440, -200, 0.01, 200, 440]
    z = [-400, 0.01, 400]
    
    # a = input('Enter increments in mm:') #increments (mm)
    # a = int(a)
    
    
    for j in y:
        if j < 300:
            # print(i+a)
            y.append(j + a)
        if y==0:
            y = 0.01
    
    pos = []
    for j in y:
        for i in x:
            for k in z:
                pos.append([i, j, k])
                # print([i, j, k])
                # print(pos)
   
    
    np.savetxt("positions" + str(a) + "mm increments.csv",pos, delimiter =",",  fmt ='% s', header = 'x (mm),y (mm),z (mm)', comments='')
    
    return

# pos(120)