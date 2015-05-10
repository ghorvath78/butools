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

def CheckDMAPRepresentation (D0, D1, prec=None):
    """
    Checks if the input matrixes define a discrete time MAP.
    
    Matrices D0 and D1 must have the same size, D0 must be a 
    transient probability matrix, D1 has only non-negative
    elements, and the rowsum of D0+D1 is 1 (up to the numerical
    precision).
    
    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the DMAP to check
    D1 : matrix, shape (M,M)
        The D1 matrix of the DMAP to check
    prec : double, optional
        Numerical precision, the default value is 1e-14
    
    Returns
    -------
    r : bool 
        The result of the check
    """

    if prec==None:
        prec=butools.checkPrecision

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

def CheckDMMAPRepresentation (D,prec=None):
    """
    Checks if the input matrixes define a discrete time MMAP.
    
    All matrices D0...DK must have the same size, D0 must be a 
    transient probability matrix, D1 has only non-negative 
    elements, and the rowsum of D0+D1+...+DK is 1 (up to the 
    numerical precision).
    
    Parameters
    ----------
    D : list/cell of matrices, length(K)
        The D0...DK matrices of the DMMAP to check
    
    Returns
    -------
    r : bool 
        The result of the check
    """

    if prec==None:
        prec=butools.checkPrecision

    if np.min(np.hstack(D)) < -prec:
        if butools.verbose:
            print ("CheckDMMAPRepresentation: Some of the matrices D1 ... DM have negative elements!")
        return False
    return CheckDMAPRepresentation(D[0],SumMatrixList(D[1:]),prec)

def CheckDRAPRepresentation (D0, D1, prec=None):
    """
    Checks if the input matrixes define a discrete time RAP.
    
    Matrices H0 and H1 must have the same size, the dominant
    eigenvalue of H0 is real and less than 1, and the rowsum of 
    H0+H1 is 1 (up to the numerical precision).
    
    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the DRAP to check
    H1 : matrix, shape (M,M)
        The H1 matrix of the DRAP to check
    prec : double, optional
        Numerical precision, the default value is 1e-14
    
    Returns
    -------
    r : bool 
        The result of the check
    """
    
    if prec==None:
        prec=butools.checkPrecision

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

def CheckDMRAPRepresentation(H,prec=None):
    """
    Checks if the input matrixes define a discrete time MRAP.
    
    All matrices H0...HK must have the same size, the dominant
    eigenvalue of H0 is real and less than 1, and the rowsum of 
    H0+H1+...+HK is 1 (up to the numerical precision).
    
    Parameters
    ----------
    H : list/cell of matrices, length(K)
        The H0...HK matrices of the DMRAP to check
    
    Returns
    -------
    r : bool 
        The result of the check
    """

    if prec==None:
        prec=butools.checkPrecision

    return CheckDRAPRepresentation(H[0],SumMatrixList(H[1:]),prec)
