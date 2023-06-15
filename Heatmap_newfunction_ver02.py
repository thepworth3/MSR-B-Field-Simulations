# -*- coding: utf-8 -*-
"""
Created on Mon May  1 11:52:14 2023

@author: thomas
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 13:48:35 2022

@author: maedeh
"""


import numpy as np
import matplotlib.pyplot as plt
from dipole import dipole
from scipy.constants import mu_0, pi
from newfunction import *
import random
import statistics
import pandas as pd
import seaborn as sb
from position_maker import pos
import time

current_time = time.time()


inc = 10  # enter your increments in mm


pos(inc)

# reading positions from csv file which is created by the position maker script
positions = pd.read_csv("positions" + str(inc) + "mm increments.csv")

x = np.array(positions.iloc[:, 0]/1000)  # Converting position to meters
# .iloc = integer location, selects rows and columns from data frame with integer based indexing (naming)
y = np.array(positions.iloc[:, 1]/1000)
z = np.array(positions.iloc[:, 2]/1000)
print(len(x))

'''
cell
'''


# positions, sets mapp grid around cell to be 1m x 1m x 1m.
xroi, yroi, zroi = np.mgrid[-.3: .3: 61j, -.3: .3: 61j, -.3: .3: 61j]

print(xroi.size)

# Cell dimensions, the nEDM cell is a cylinder whichs stands up vertically inside of the MSR
rcell = 0.18  # m, cell radius
hcell = 0.15  # m, cell height
dcell = 0.10  # m, bottom to top distance of cells


mask = (abs(zroi) >= dcell/2) & (abs(zroi) <= dcell /
                                 2 + hcell) & (xroi ** 2 + yroi ** 2 < rcell ** 2)

xroi, yroi, zroi = xroi[mask], yroi[mask], zroi[mask]

print(xroi.size)

l1 = int(input('Enter your desire order:'))
l2 = int(input('Enter your desire order to fit:'))


Hx1, Hy1, Hz1 = CreateBListToOrder(l1)     # b = mu*H
Hx2, Hy2, Hz2 = CreateBListToOrder(l2)

'''
Munich glm B fields
'''

MunicGlm = glmCreator(l1)     # called from newfunction code



k = x.size
#print("Generating the Munich glm took", time.time()-glm_time)
residual_Bx_total = []  # empty lsit for residuals of B
residual_By_total = []
residual_Bz_total = []

randomlistx = []  # empyt list for positions x, y, z
randomlisty = []
randomlistz = []

randomlistBx = []  # empty list for B field values
randomlistBy = []
randomlistBz = []

# standard deviation of the position arguements are (# poserror, # Berror)
stdevx = np.zeros((5, 7))
stdevy = np.zeros((5, 7))
stdevz = np.zeros((5, 7))

stdevB = np.zeros((5, 7))

stdevxx = []  # standard devaiton for each component of the radius rr in spherical coords
stdevyy = []
stdevzz = []

# different positions errors that we want to simualte for
poserror = (0, 0.005, 0.01, 0.015, 0.02)
# different field errors that we want to simualte for
Berror = (.1e-12, 0.5e-12, 1e-12, 5e-12, 10e-12, 50e-12, 100e-12)

# u and v are indices for each value of pos error and B error that we have so we can label them
u = (1, 2, 3, 4, 5)
v = (1, 2, 3, 4, 5, 6, 7)
a = 0
b = 0

rangej = 500 # how many sets of random numbers we want to generate

for countpe, pe in enumerate(poserror):               #
    for countbe, Be in enumerate(Berror):
        residual_Bx_total = [0]
        residual_By_total = [0]
        residual_Bz_total = [0]
        residual_B_total = [0]

        for j in range(rangej):

            # attach a randomly generated error to the x coord
            randomlistx = np.array([random.uniform(-pe, pe) for i in range(k)])
            Xwitherror = x + randomlistx

            # attach a randomly generated error to the y coord
            randomlisty = np.array([random.uniform(-pe, pe) for i in range(k)])
            Ywitherror = y + randomlisty

            # attach a randomly generated error to the z coord
            randomlistz = np.array([random.uniform(-pe, pe) for i in range(k)])
            Zwitherror = z + randomlistz

            BX = Bx(Hx1, x, y, z, MunicGlm)
            # for every set of values of x,y,z generate the Bx error value and append to a list
            randomlistBx = np.array(
                [np.random.normal(loc=0.0, scale=Be, size=None) for i in range(k)])
            BXX = BX + randomlistBx  # add the random Bx error to the calculated Bx using Munich Glm

            BY = By(Hy1, x, y, z, MunicGlm)
            randomlistBy = np.array(
                [np.random.normal(loc=0.0, scale=Be, size=None) for i in range(k)])
            BYY = BY + randomlistBy

            BZ = Bz(Hz1, x, y, z, MunicGlm)
            randomlistBz = np.array(
                [np.random.normal(loc=0.0, scale=Be, size=None) for i in range(k)])
            BZZ = BZ + randomlistBz

            Mfitresult = glm_fit(Hx2, Hy2, Hz2, l2, Xwitherror,
                                  Ywitherror, Zwitherror, BXX, BYY, BZZ)
            # Mresiduals = Mfitresult[1]     # fit the result to get the new simulated glm for our MSR
            newglm = Mfitresult[0]

            # print(newglm)

            # now lets cacluate the Bx , By, Bz over the xroi, yroi, zroi
            
            Bxroi = Bx(Hx2, xroi, yroi, zroi, MunicGlm)
            Byroi = By(Hy2, xroi, yroi, zroi, MunicGlm)
            Bzroi = Bz(Hz2, xroi, yroi, zroi, MunicGlm)
            Broi = np.sqrt(Bxroi ** 2 + Byroi ** 2 + Bzroi ** 2)

            # # now lets cacluate the Bx , By, Bz over the xroi, yroi, zroi with the new glms

            Bxroi_err = Bx(Hx2, xroi, yroi, zroi, newglm)
            Byroi_err = By(Hy2, xroi, yroi, zroi, newglm)
            Bzroi_err = Bz(Hz2, xroi, yroi, zroi, newglm)
            Broi_err = np.sqrt(Bxroi_err ** 2 + Byroi_err **2 + Bzroi_err ** 2)

            #now we calculate the difference (residual over the ROI)
            
            residual_Bx = np.squeeze([Bxroi - Bxroi_err])
            residual_By = np.squeeze([Byroi - Byroi_err])
            residual_Bz = np.squeeze([Bzroi - Bzroi_err])
            
            residual_B = np.squeeze([Broi - Broi_err])
         
            
            residual_Bx_total.extend(residual_Bx)
            residual_By_total.extend(residual_By)
            residual_Bz_total.extend(residual_Bz)
            
            residual_B_total.extend(residual_B)
            print(j)
            
        stdevx[countpe][countbe] = (
             np.std(residual_Bx_total))*(10e12)
        stdevy[countpe][countbe] = (
             np.std(residual_By_total))*(10e12)
        stdevz[countpe][countbe] = (
             np.std(residual_Bz_total))*(10e12)
        stdevB[countpe][countbe] = (np.std(residual_B_total))*(10e12)
        print("B err number",countbe+1)
    print("Pos err number", countpe+1 )
nBerror = (.1, 0.5, 1, 5, 10, 50, 100)


# # Bx heatmap

# # dataframe is a 2D labeled structure with columns of different types (each column is different order)
dfx = pd.DataFrame(data=stdevx, index=[*poserror], columns=[*nBerror])
dfx = dfx.rename_axis("position error (m)", axis="index")
dfx = dfx.rename_axis("B error (pT)", axis="columns")


fig, ax = plt.subplots()
sb.heatmap(dfx, annot=True, cmap='plasma', fmt=".2f", cbar_kws={
            'label': 'Std. Dev. of Residuals in Bx (pT)'})
ax.invert_yaxis()

plt.savefig("NewFunction_pebe_" + str(rangej) + "_" + str(inc) +
             "mm - Bx_cell.png", bbox_inches='tight', dpi=320)


# # By heatmap

dfy = pd.DataFrame(data=stdevy, index=[*poserror], columns=[*nBerror])
# dataframe is a 2D labeled structure with columns of different types (each column is different order)
dfy = dfy.rename_axis("position error (m)", axis="index")
dfy = dfy.rename_axis("B error (pT)", axis="columns")

fig2, ax2 = plt.subplots()
sb.heatmap(dfy, annot=True, cmap='plasma', fmt=".2f", cbar_kws={
            'label': 'Std. Dev. of Residuals in By (pT)'})
ax2.invert_yaxis()

plt.savefig("NewFunction_pebe_" + str(rangej) + "_" + str(inc) +
             "mm - By_cell.png", bbox_inches='tight', dpi=320)


# # Bz haetmap

# # dataframe is a 2D labeled structure with columns of different types (each column is different order)
dfz = pd.DataFrame(data=stdevz, index=[*poserror], columns=[*nBerror])
dfz = dfz.rename_axis("position error (m)", axis="index")
dfz = dfz.rename_axis("B error (pT)", axis="columns")

fig3, ax3 = plt.subplots()
sb.heatmap(dfz, annot=True, cmap='plasma', fmt=".2f", cbar_kws={
            'label': 'Std. Dev. of Residuals in Bz (pT)'})
ax3.invert_yaxis()

plt.savefig("NewFunction_pebe_" + str(rangej) + "_" + str(inc) +
             "mm - Bz_cell.png", bbox_inches='tight', dpi=320)


# # B tot heatmap

dfB = pd.DataFrame(data=stdevB, index=[*poserror], columns=[*nBerror])
dfB = dfB.rename_axis("position error (m)", axis="index")
dfB = dfB.rename_axis("B error (pT)", axis="columns")

fig1, ax1 = plt.subplots()
sb.heatmap(dfB, annot=True, cmap='plasma', fmt=".2f", cbar_kws={
           'label': 'Std. Dev. of Residuals in B (pT)'})
ax1.invert_yaxis()

plt.savefig("NewFunction_pebe_" + str(rangej) + "_" + str(inc) +
             "mm - Btot_cell.png", bbox_inches='tight', dpi=320)


# # Histogram

fig4, ax4 = plt.subplots()

n, bins, patches = ax4.hist(BZ, color='purple')
ax4.set_title(r'Histogram $B_{z}$ over the cell & door holes')
ax4.set_xlabel('$\Delta B_z$ (T)')
ax4.set_ylabel('# of entries')
ax4.tick_params(axis='both')
fig4.savefig(r'NEW Histogram $B_{z}$ over the cell & door holes.png', dpi=320)

plt.show()

#Bz vs. z

fig5, ax5 = plt.subplots()

ax5.plot(z, BZ, 'r.', label='$B_z$')
#ax1.plot(zroi[mask], BzM_new - BzM_fit, '.', label='differences $B_z$')
#ax1.plot(zroi[mask], BzM_fit, 'b.', label='fit $B_z$')
ax5.set_xlabel(r'z (m)')
ax5.set_ylabel(r'$B_z (T)$')
ax5.set_title(r'$B_z$ fields (Munich $G_{lm}$) & increments = ' + str(inc) +
               'mm over the cell', y=1.04, fontsize=10)
ax5.legend()
fig5.savefig(r'Hole pos $B_z$ fields Munich glm & increments = '
              + str(inc) + 'mm.png', dpi=320)


fig1 = plt.figure()

ax1 = fig1.add_subplot(111)
n, bins, patches = ax1.hist(residual_Bx_total, bins = 500, color = 'navy')
ax1.set_title('1')
ax1.set_xlabel('residual Bx (T)')
ax1.set_ylabel('# of entries')

plt.savefig("randx1WE.png",bbox_inches='tight', dpi=320)


fig2 = plt.figure()

ax2 = fig2.add_subplot(111)
n, bins, patches = ax2.hist(residual_By_total, bins = 500, color = 'green')
ax2.set_title('1')
ax2.set_xlabel('residual By (T)')
ax2.set_ylabel('# of entries')

plt.savefig("randy1WE.png",bbox_inches='tight', dpi=320)


fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
n, bins, patches = ax3.hist(residual_Bz_total, bins = 500, color = 'orange')
ax3.set_title('1')
ax3.set_xlabel('residual Bz (T)')
ax3.set_ylabel('# of entries')

plt.savefig("randz1WE.png",bbox_inches='tight', dpi=320)


fig4 = plt.figure()

ax4 = fig4.add_subplot(111)
n, bins, patches = ax4.hist(residual_B_total, bins = 500, color = 'purple')
ax4.set_title('1')
ax4.set_xlabel('residual B (T)')
ax4.set_ylabel('# of entries')

plt.savefig("rand1WE.png",bbox_inches='tight', dpi=320)



process_time = time.time()-current_time

print("The simulation took",  process_time, "seconds")

