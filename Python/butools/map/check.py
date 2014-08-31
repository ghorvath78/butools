# -*- coding: utf-8 -*-
"""
Created on Sun Aug 24 15:31:04 2014

@author: gabor
"""

import butools
from butools.mc import CheckGenerator, CheckProbVector
from butools.utils import SumMatrixList
import numpy as np
import scipy.linalg as la

def CheckMAPRepresentation (D0, D1, prec=1e-14):

    if not CheckGenerator(D0,True):
        return False

    if D0.shape!=D0.shape:
        if butools.verbose:
            print ("CheckMAPRepresentation: D0 and D1 have different sizes!")
        return False

    if np.min(D1)<-prec:
        if butools.verbose:
            print ("CheckMAPRepresentation: D1 has negative element!")
        return False

    if np.any(np.abs(np.sum(D0+D1,1))>prec):
        if butools.verbose:
            print ("CheckMAPRepresentation: The rowsum of D0+D1 is not 0!")
        return False

    return True

def CheckMMAPRepresentation (H,prec=1e-14):

    if np.min(np.hstack(H[1:])) < -prec:
        if butools.verbose:
            print ("CheckMMAPRepresentation: Some of the matrices H1 ... HM have a negative element!")
        return False
    return CheckMAPRepresentation(H[0],SumMatrixList(H[1:]),prec)

def CheckRAPRepresentation (D0, D1, prec=1e-14):

    if D0.shape[0]!=D0.shape[1]:
        if butools.verbose:
            print ("CheckRAPRepresentation: D0 is not a quadratic matrix!")
        return False

    if D1.shape[0]!=D1.shape[1]:
        if butools.verbose:
            print ("CheckRAPRepresentation: D1 is not a quadratic matrix!")
        return False

    if D0.shape!=D1.shape:
        if butools.verbose:
            print ("CheckRAPRepresentation: D0 and D1 have different sizes!")
        return False

    if np.max(np.abs((D0+D1).dot(np.ones((D0.shape[0],1))))) > prec:
        if butools.verbose:
            print ("CheckRAPRepresentation: A rowsum of D0+D1 is not 0!")
        return False

    ev=la.eigvals(D0)
    if np.max(np.real(ev))>=-prec:
        if butools.verbose:
            print ("heckRAPRepresentation: there is an eigenvalue of D0 with non-negative real part")
        return False

    ceig=np.array(ev)
    reig=np.array(ev)
    for i in range(len(ev)):
        if np.isreal(ev[i]):
            ceig[i]=-np.Inf
            reig[i]=np.real(ev[i])
        else:
            ceig[i]=np.real(ev[i])
            reig[i]=-np.Inf;

    if np.max(reig) < np.max(ceig):
        if butools.verbose:
            print ("CheckRAPRepresentation: The dominant eigenvalue of D0 is not real!")
        return False

    if np.max(reig)==np.max(ceig):
        if butools.verbose:
            print("CheckRAPRepresentation: The dominant and a complex eigenvalue of D0 has the same real part!")

    return True

def CheckMRAPRepresentation(H,prec=1e-14):

    return CheckRAPRepresentation(H[0],SumMatrixList(H[1:]),prec)
