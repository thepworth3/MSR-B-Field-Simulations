# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 14:05:23 2023

@author: thoma
"""
import numpy as np
import matplotlib.pyplot as plt
from dipole import dipole
from scipy.constants import mu_0, pi
from newfunction_v1 import *
import random
import statistics
import pandas as pd
import seaborn as sb
from position_maker_dipole_sim import pos
import time

inc = 10 #[10, 50, 100, 200, 500] # enter your increments in mm

center_of_MSR = False

  
pos(inc)
positions = pd.read_csv("positions" + str(inc) + "mm increments.csv")
x = np.array([positions.iloc[:, 0]/1000])  # Converting position to meters
y = np.array([positions.iloc[:, 1]/1000])
z = np.array([positions.iloc[:, 2]/1000])


x = x.T.reshape(x.size)
y = y.T.reshape(x.size) #Gets the dimensions of the positiona arrays to be proper, allowing the T matrix to compute properly
z = z.T.reshape(x.size)

# '''
# cell
# '''

# # positions, sets mapp grid around cell to be 1m x 1m x 1m.
# xroi, yroi, zroi = np.mgrid[-.3: .3: 61j, -.3: .3: 61j, -.3: .3: 61j]


# # Cell dimensions, the nEDM cell is a cylinder whichs stands up vertically inside of the MSR
# rcell = 0.18  # m, cell radius
# hcell = 0.15  # m, cell height
# dcell = 0.10  # m, bottom to top distance of cells


# mask = (abs(zroi) >= dcell/2) & (abs(zroi) <= dcell /
#                                   2 + hcell) & (xroi ** 2 + yroi ** 2 < rcell ** 2)

# xroi, yroi, zroi = xroi[mask], yroi[mask], zroi[mask]


l1 = int(input('Enter your desire order:'))
l2 = int(input('Enter your desire order to fit:'))


Hx1, Hy1, Hz1 = CreateBListToOrder(l1)     # b = mu*H
Hx2, Hy2, Hz2 = CreateBListToOrder(l2)

'''
Munich glm B fields
'''

MunicGlm = glmCreator(l1)     # called from newfunction code


xd = 0
yd = 1.2                      #manually defining the position of the dipole
zd = 0

mx = 0
my = 0.01                               #manually defining the magnetic moment of dipole 
mz = 0

# d = dipole()
d = dipole()
d.set(xd, yd, zd, mx, my, mz)
# d = dipole(xd, yd, zd, mx, my, mz)

Bxd = d.bx(x, y, z)
Byd = d.by(x, y, z)           #B calc with dipole code
Bzd = d.bz(x, y, z)
Btd = np.sqrt(Bxd**2 + Byd**2 + Bzd**2)

# print("Shape of Hx1",Hx1.shape)
# print("Shape of MunicGlm", MunicGlm.shape)
# print("Shape of x", x.shape)
# print("Type of x", x.dtype)



Mfitresult = glm_fit(Hx2, Hy2, Hz2, l2, x, y, z, Bxd, Byd, Bzd)

# Mresiduals = Mfitresult[1]     # fit the result to get the new simulated glm for our MSR

newglm = Mfitresult[0]
BXG = Bx(Hx1, x, y, z, newglm)
BYG = By(Hy1, x, y, z, newglm)
BZG = Bz(Hz1, x, y, z, newglm)
BTG = B(Hx1,Hy1,Hz1, x, y, z, newglm)

residualx = BXG - Bxd
residualy = BYG - Byd
residualz = BZG - Bzd
residualt = BTG - Btd


stdevx = np.std(residualx)
stdevy = np.std(residualy)
stdevz = np.std(residualz)
stdevt = np.std(residualt)
print(stdevt)
# print("new glms are", newglm)


# filename = "Simulations_data_order" + str(l1) + "_inc" + str(inc) 
# df = pd.DataFrame()
# pd.DataFrame.to_csv(filename, sep=',', index = False, encoding='utf-8')




#Now we can calculate the field with the new glms so we can compare

with plt.style.context("bmh"):
    # plt.plot(x, residualx, ".g", label = "Residual Bx")
    # plt.plot(x, residualy, ".r", label = "Residual By")
    # plt.plot(x, residualz, ".b", label = "Residual Bz")
    plt.plot(x, residualt, ".y", label = "Residual Btot")
    plt.xlabel("x position (m)")
    plt.ylabel("Residual B field (T)")
    plt.legend()
    
    plt.show()
    
    # plt.plot(y, residualx, ".g", label = "Residual Bx")
    # plt.plot(y, residualy, ".r", label = "Residual By")
    # plt.plot(y, residualz, ".b", label = "Residual Bz")
    plt.plot(y, residualt, ".y", label = "Residual Btot")
    plt.xlabel("y position (m)")
    plt.ylabel("Residual B field (T)")
    plt.legend()
    
    plt.show()
    
    # plt.plot(z, residualx, ".g", label = "Residual Bx")
    # plt.plot(z, residualy, ".r", label = "Residual By")
    # plt.plot(z, residualz, ".b", label = "Residual Bz")
    plt.plot(z, residualt, ".y", label = "Residual Btot")
    plt.xlabel("z position (m)")
    plt.ylabel("Residual B field (T)")
    plt.legend()
    
    plt.show()
    
    plt.hist(residualt, bins = 50)
    plt.xlabel("Residual Magnetic field (T) ")
    plt.ylabel("Number of entries ")
    plt.show()

plt.plot(y,Bzd, "go",label ="dipole")
plt.plot(y,BZG, "r.",label ="glm fit")
plt.xlabel("y position (m)")
plt.ylabel("Bz (T)")
plt.legend()
plt.show()

plt.plot(x,Bzd, "g.",label ="dipole")
plt.plot(x,BZG, "r.",label ="glm fit")
plt.xlabel("x position (m)")
plt.ylabel("Bz (T)")
plt.legend()
plt.show()

plt.plot(z,Bzd, "g.",label ="dipole")
plt.plot(z,BZG, "r.",label ="glm fit")
plt.xlabel("z position (m)")
plt.ylabel("Bz (T)")
plt.legend()
plt.show()






# plt.plot(x, Bxd, "ro", label = "Dipole field")
# plt.plot(x, BXG, "g.", label = "glm field")
# plt.plot(x, percent_diff, "r.")
# plt.xlabel("position")
# plt.ylabel("B field")
# plt.legend()
# plt.title("B field comparison for order =", l1)
# plt.show()








'''

# Plotting  3D
fig31 = plt.figure()
ax31 = plt.axes(projection = "3d")

plt.title("Bx over the xy plane ")
ax31.scatter3D(xroi,yroi,Bxd)
ax31.set_xlabel("x-axis (m)")
ax31.set_ylabel("y-axis (m)")
ax31.set_zlabel("Bx from dipole (T)") # Need to fix this scaling, ask Maedeh is she knows how
plt.autoscale(False)
plt.show()

fig32 = plt.figure()
ax32 = plt.axes(projection = "3d")

plt.title("By over the xy plane")
ax32.scatter3D(x,y,Byd)
ax32.set_xlabel("x-axis (m)")
ax32.set_ylabel("y-axis (m)")
ax32.set_zlabel("By from dipole (T)") # Need to fix this scaling, ask Maedeh is she knows how
plt.autoscale(False)
plt.show()

fig33 = plt.figure()
ax33 = plt.axes(projection = "3d")

plt.title("Bz over the xy plane")
ax33.scatter3D(x,y,Bzd)
ax33.set_xlabel("x-axis (m)")
ax33.set_ylabel("y-axis (m)")
ax33.set_zlabel("Bz from dipole (T)") # Need to fix this scaling, ask Maedeh is she knows how
plt.autoscale(False)
plt.show()


fig34 = plt.figure()
ax34= plt.axes(projection = "3d")

plt.title("Bx over the xz plane")
ax34.scatter3D(x,z,Bxd)
ax34.set_xlabel("x-axis (m)")
ax34.set_ylabel("z-axis (m)")
ax34.set_zlabel("Bx from dipole (T)") # Need to fix this scaling, ask Maedeh is she knows how
plt.autoscale(False)
plt.show()

fig35 = plt.figure()
ax35 = plt.axes(projection = "3d")

plt.title("By over the xz plane")
ax35.scatter3D(x,z,Byd)
ax35.set_xlabel("x-axis (m)")
ax35.set_ylabel("z-axis (m)")
ax35.set_zlabel("By from dipole (T)") # Need to fix this scaling, ask Maedeh is she knows how
plt.autoscale(False)
plt.show()

fig36 = plt.figure()
ax36 = plt.axes(projection = "3d")

plt.title("Bz over the xz plane")
ax36.scatter3D(x,z,Bzd)
ax36.set_xlabel("x-axis (m)")
ax36.set_ylabel("z-axis (m)")
ax36.set_zlabel("Bz from dipole (T)") # Need to fix this scaling, ask Maedeh is she knows how
plt.autoscale(False)
plt.show()

fig37 = plt.figure()
ax37 = plt.axes(projection = "3d")

plt.title("Bx over the yz plane")
ax37.scatter3D(y,z,Bxd)
ax37.set_xlabel("y-axis (m)")
ax37.set_ylabel("z-axis (m)")
ax37.set_zlabel("Bx from dipole (T)") # Need to fix this scaling, ask Maedeh is she knows how
plt.autoscale(False)
plt.show()


fig36 = plt.figure()
ax36 = plt.axes(projection = "3d")

plt.title("By over the yz plane")
ax36.scatter3D(y,z,Byd)
ax36.set_xlabel("y-axis (m)")
ax36.set_ylabel("z-axis (m)")
ax36.set_zlabel("By from dipole (T)") # Need to fix this scaling, ask Maedeh is she knows how
plt.autoscale(False)
plt.show()

fig36 = plt.figure()
ax36 = plt.axes(projection = "3d")

plt.title("Bz over the yz plane")
ax36.scatter3D(y,z,Bzd)
ax36.set_xlabel("y-axis (m)")
ax36.set_ylabel("z-axis (m)")
ax36.set_zlabel("Bz from dipole (T)") # Need to fix this scaling, ask Maedeh is she knows how
plt.autoscale(False)
plt.show()


# to get any thing we can fit, need to adjust how the positions work and just do one hole at a time
#Plotting 2D    
fig21 = plt.figure()
ax21 = plt.axes()

plt.title("Bx vs. x position")
plt.plot(x,Bxd,"ro")
plt.xlabel('x axis (m)')
plt.ylabel('Bx due to dipole (T)')



fig22 = plt.figure()
ax22 = plt.axes()

plt.title("By vs. x position")
plt.plot(x,Byd,"ro")
plt.xlabel('x axis (m)')
plt.ylabel('By due to dipole (T)')

fig23 = plt.figure()
ax23 = plt.axes()

plt.title("Bz vs. x position")
plt.plot(x,Bzd,"ro")
plt.xlabel('x axis (m)')
plt.ylabel('Bz due to dipole (T)')


fig24 = plt.figure()
ax24 = plt.axes()

plt.title("Bx vs. y position")
plt.plot(y,Bxd,"ro")
plt.xlabel('y axis (m)')
plt.ylabel('Bx due to dipole (T)')


fig25 = plt.figure()
ax25 = plt.axes()

plt.title("By vs. y position")
plt.plot(y,Byd,"ro")
# plt.plot(x+0.5, -1e-10/x**3)
plt.xlabel('y axis (m)')
plt.ylabel('Bz due to dipole (T)')

fig26 = plt.figure()
ax26 = plt.axes()

plt.title("Bz vs. y position")
plt.plot(y,Bzd,"ro")
plt.xlabel('y axis (m)')
plt.ylabel('Bz due to dipole (T)')

fig27 = plt.figure()
ax27 = plt.axes()

plt.title("Bx vs. z position")
plt.plot(z,Bxd, "co")
plt.xlabel('z axis (m)')
plt.ylabel('Bx due to dipole (T)')

fig28 = plt.figure()
ax28 = plt.axes()

plt.title("By vs. z position")
plt.plot(z,Byd, "co")
plt.xlabel('y axis (m)')
plt.ylabel('By due to dipole (T)')

fig27 = plt.figure()
ax27 = plt.axes()

plt.title("Bz vs. z position")
plt.plot(z,Bzd, "co")
plt.xlabel('z axis (m)')
plt.ylabel('Bz due to dipole (T)')

fig28 = plt.figure()
ax28 = plt.axes()

plt.title("Btot vs. x")
plt.plot(x,Btd, "go")
plt.xlabel('x axis (m)')
plt.ylabel('Btot due to dipole (T)')

fig29 = plt.figure()
ax29 = plt.axes()

plt.title("Btot vs. y")
plt.plot(y,Btd, "go")
plt.xlabel('y axis (m)')
plt.ylabel('Btot due to dipole (T)')

fig291 = plt.figure()
ax291 = plt.axes()

plt.title("Btot vs. z")
plt.plot(z,Btd, "go")
plt.xlabel('z axis (m)')
plt.ylabel('Btot due to dipole (T)')
'''