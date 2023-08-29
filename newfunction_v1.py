#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 15:01:49 2022

@author: maedeh
v1 by Thomas Hepworth July 2023. Removed reshaping of the B functions
"""
from PislibHelperFunctions import *
import numpy as np
from scipy.linalg import lstsq


# c1 = 0.7957747150262e-6        #to change Glm in pT/cm^l to A/m^(l+1)
c1 = 1e-12
G1 = c1 * 30
G2 = c1 * 6e-1
G3 = c1 * 6e-3
G4 = c1 * 9.5e-5
G5 = c1 * 9e-7


# creating T matrix

def TxValue(Hx, x, y, z):

    Tx = []

    for t in Hx:
        Tx.append(t(x, y, z))
        if np.isscalar(Tx[-1]):
            Tx[-1] = np.full_like(x, Tx[-1])
    Tx = np.array(Tx).T

    return(Tx)


def TyValue(Hy, x, y, z):

    Ty = []

    for t in Hy:
        Ty.append(t(x, y, z))
        if np.isscalar(Ty[-1]):
            Ty[-1] = np.full_like(y, Ty[-1])
    Ty = np.array(Ty).T

    return(Ty)


def TzValue(Hz, x, y, z):

    Tz = []

    for t in Hz:
        Tz.append(t(x, y, z))
        if np.isscalar(Tz[-1]):
            Tz[-1] = np.full_like(z, Tz[-1])
    Tz = np.array(Tz).T

    return(Tz)


def HValue(Hx, Hy, Hz, x, y, z):

    Tx = TxValue(Hx, x, y, z)
    Ty = TyValue(Hy, x, y, z)
    Tz = TzValue(Hz, x, y, z)

    T = np.concatenate((Tx, Ty, Tz))

    return T


def glmCreator(l):

    G = [0, G1, G2, G3, G4, G5]

    gt7 = []

    for l in range(l+1):
        glm = G[l]
        k = 2*(l+1)+1
        gt7.extend([glm]*k)
        l = l+1

    gt7 = np.asarray(gt7)
    # .reshape(len(gt7),1)
    return (gt7)


# Bx
def Bx(Hx, x, y, z, glms):
    Bx = TxValue(Hx, x, y, z).dot(glms)
    # Bx = np.asarray(Bx).reshape(len(x))
    return Bx


# By
def By(Hy, x, y, z, glms):
    By = TyValue(Hy, x, y, z).dot(glms)
    # By = np.asarray(By).reshape(len(x))
    return By


# Bz
def Bz(Hz, x, y, z, glms):
    Bz = TzValue(Hz, x, y, z).dot(glms)
    # Bz = np.asarray(Bz).reshape(len(x))
    return Bz


# True B field
def B(Hx, Hy, Hz, x, y, z, glms):
    Bx = TxValue(Hx, x, y, z).dot(glms)
    Bx = np.asarray(Bx).reshape(len(x))
    By = TyValue(Hy, x, y, z).dot(glms)
    By = np.asarray(By).reshape(len(x))
    Bz = TzValue(Hz, x, y, z).dot(glms)
    Bz = np.asarray(Bz).reshape(len(x))
    B = np.sqrt(Bx**2 + By**2 + Bz**2)
    return B


# fit glm
def glm_fit(Hx, Hy, Hz, l, x, y, z, Bx, By, Bz):
    T = HValue(Hx, Hy, Hz, x, y, z)
    B = np.concatenate((Bx, By, Bz))
    print("Shape of T is ", T.shape)
    print("Shape of B is ", B.shape)
    gfit = lstsq(T, B)

    return gfit
