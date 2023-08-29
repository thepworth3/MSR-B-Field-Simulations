# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 21:34:14 2023

@author: thoma
"""

import matplotlib.pyplot as plt


sigma = [2.0249e-11, 2.1791e-11, 2.3292e-11, 2.6769e-11,]# 3.2496e-11, 3.8344e-11,1.8837e-10,4.3344e-10, 1.0340e-09]
inc = [1,3,5,10]#,20,30,40,45,50]

with plt.style.context("seaborn-v0_8-deep"):
    plt.plot(inc,sigma, linestyle="-", marker="o",color="g")
    plt.xlabel("Measurement increment (cm)")
    plt.ylabel("Standard Deviation of Residuals (T)")

