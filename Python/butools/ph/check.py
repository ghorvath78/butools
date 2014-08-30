# -*- coding: utf-8 -*-
"""
Created on Sun Aug 24 15:31:04 2014

@author: gabor
"""

import butools
from butools.mc import CheckGenerator, CheckProbVector
import numpy as np
import scipy.linalg as la

def CheckPHRepresentation (alpha, A, prec=1e-14):

    if len(alpha.shape)<2:
        if butools.verbose:
            print("CheckPHRepresentation: Initial vector must be a matrix of size (1,N)!")
        return False
    
    if alpha.shape[1]!=A.shape[0]:
        if butools.verbose:
            print("CheckPHRepresentation: The vector and the matrix have different sizes!")
        return False

    if not CheckGenerator(A,True,prec) or not CheckProbVector(alpha,True,prec):
        return False

    return True

def CheckMERepresentation (alpha, A, prec=1e-14):

    if len(alpha.shape)<2:
        if butools.verbose:
            print("CheckMERepresentation: Initial vector must be a matrix of size (1,N)!")
        return False

    if A.shape[0]!=A.shape[1]:
        if butools.verbose:
            print("CheckMERepresentation: The matrix is not a square matrix!")
        return False

    if alpha.shape[1]!=A.shape[0]:
        if butools.verbose:
            print("CheckMERepresentation: The vector and the matrix have different sizes!")
        return False

    if np.sum(alpha)<-prec*alpha.size or np.sum(alpha)>1.0+prec*alpha.size:
        if butools.verbose:
            print("CheckMERepresentation: The sum of the vector elements is less than zero or greater than one!")
        return False

    ev = la.eigvals(A)

    if np.max(np.real(ev))>=prec:
        if butools.verbose:
            print("CheckMERepresentation: There is an eigenvalue of the matrix with non-negative real part!")
        return False
    
    ix = np.argsort(np.abs(np.real(ev)))
    maxev = ev[ix[0]]

    if not np.isreal(maxev):
        if butools.verbose:
            print("CheckMERepresentation: The dominant eigenvalue of the matrix is not real!")
        return False

    if np.sum(np.abs(ev)==abs(maxev)) > 1 and butools.verbose:
        print ("CheckMERepresentation warning: There are more than one eigenvalue with the same absolute value as the largest eigenvalue!")

    return True
    
