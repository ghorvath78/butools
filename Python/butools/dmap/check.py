# -*- coding: utf-8 -*-
"""
Created on Sun Aug 24 15:31:04 2014

@author: gabor
"""

import butools
from butools.mc import CheckProbMatrix
from butools.utils import SumMatrixList
import numpy as np
import scipy.linalg as la

def CheckDMAPRepresentation (D0, D1, prec=1e-14):

    if not CheckProbMatrix(D0,True, prec):
        if butools.verbose:
            print ("CheckDMAPRepresentation: D0 is not a transient probability matrix!")
        return False

    if D0.shape!=D1.shape:
        if butools.verbose:
            print ("CheckDMAPRepresentation: D0 and D1 have different sizes!")
        return False

    if np.min(D1)<-prec or np.min(D0)<-prec:
        if butools.verbose:
            print ("CheckDMAPRepresentation: One of the matrices has negative element!")
        return False

    if np.any(np.abs(np.sum(D0+D1,1)-1)>prec):
        if butools.verbose:
            print ("CheckDMAPRepresentation: The rowsum of D0+D1 is not 1!")
        return False

    return True

def CheckDMMAPRepresentation (D,prec=1e-14):

    if np.min(np.hstack(D)) < -prec:
        if butools.verbose:
            print ("CheckDMMAPRepresentation: Some of the matrices D1 ... DM have negative elements!")
        return False
    return CheckDMAPRepresentation(D[0],SumMatrixList(D[1:]),prec)

def CheckDRAPRepresentation (D0, D1, prec=1e-14):

    if D0.shape[0]!=D0.shape[1]:
        if butools.verbose:
            print ("CheckDRAPRepresentation: D0 is not a quadratic matrix!")
        return False

    if D1.shape[0]!=D1.shape[1]:
        if butools.verbose:
            print ("CheckDRAPRepresentation: D1 is not a quadratic matrix!")
        return False

    if D0.shape!=D1.shape:
        if butools.verbose:
            print ("CheckDRAPRepresentation: D0 and D1 have different sizes!")
        return False

    if np.any(np.abs(np.sum(D0+D1,1).A.flatten()-1.0) > prec):
        if butools.verbose:
            print ("CheckDRAPRepresentation: A rowsum of D0+D1 is not 1!")
        return False

    ev = la.eigvals(D0)
    ix = np.argsort(-np.abs(np.real(ev)))
    maxev = ev[ix[0]]

    if not np.isreal(maxev):
        if butools.verbose:
            print("CheckDRAPRepresentation: The largest eigenvalue of matrix D0 is complex!")
        return False

    if maxev>1.0+prec:
        if butools.verbose:
            print("CheckDRAPRepresentation: The largest eigenvalue of matrix D0 is greater than 1!")
        return False       

    if np.sum(np.abs(ev)==abs(maxev)) > 1 and butools.verbose:
        print ("CheckDRAPRepresentation warning: There are more than one eigenvalue with the same absolute value as the largest eigenvalue!")

    return True

def CheckDMRAPRepresentation(H,prec=1e-14):

    return CheckDRAPRepresentation(H[0],SumMatrixList(H[1:]),prec)
