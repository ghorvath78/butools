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

def CheckMAPRepresentation (D0, D1, prec=None):
    """
    Checks if the input matrixes define a continuous time MAP.
    
    Matrices D0 and D1 must have the same size, D0 must be a 
    transient generator matrix, D1 has only non-negative 
    elements, and the rowsum of D0+D1 is 0 (up to the numerical
    precision).
    
    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the MAP to check
    D1 : matrix, shape (M,M)
        The D1 matrix of the MAP to check
    prec : double, optional
        Numerical precision, the default value is 1e-14
    
    Returns
    -------
    r : bool 
        The result of the check
    """

    if prec==None:
        prec=butools.checkPrecision

    if not CheckGenerator(D0,True):
        return False

    if D0.shape!=D1.shape:
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

def CheckMMAPRepresentation (H,prec=None):
    """
    Checks if the input matrixes define a continuous time MMAP.
    
    All matrices D0...DK must have the same size, D0 must be a 
    transient generator matrix, D1 has only non-negative 
    elements, and the rowsum of D0+D1+...+DK is 0 (up to the 
    numerical precision).
    
    Parameters
    ----------
    D : list/cell of matrices, length(K)
        The D0...DK matrices of the MMAP to check
    
    Returns
    -------
    r : bool 
        The result of the check
    """

    if prec==None:
        prec=butools.checkPrecision

    if np.min(np.hstack(H[1:])) < -prec:
        if butools.verbose:
            print ("CheckMMAPRepresentation: Some of the matrices H1 ... HM have a negative element!")
        return False
    return CheckMAPRepresentation(H[0],SumMatrixList(H[1:]),prec)

def CheckRAPRepresentation (D0, D1, prec=None):
    """
    Checks if the input matrixes define a continuous time RAP.
    
    Matrices H0 and H1 must have the same size, the dominant
    eigenvalue of H0 is negative and real, and the rowsum of 
    H0+H1 is 0 (up to the numerical precision).
    
    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the RAP to check
    H1 : matrix, shape (M,M)
        The H1 matrix of the RAP to check
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
            print ("CheckRAPRepresentation: there is an eigenvalue of D0 with non-negative real part")
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

def CheckMRAPRepresentation(H,prec=None):
    """
    Checks if the input matrixes define a continuous time MRAP.
    
    All matrices H0...HK must have the same size, the dominant
    eigenvalue of H0 is negative and real, and the rowsum of 
    H0+H1+...+HK is 0 (up to the numerical precision).
    
    Parameters
    ----------
    H : list/cell of matrices, length(K)
        The H0...HK matrices of the MRAP to check
    
    Returns
    -------
    r : bool 
        The result of the check
    """

    if prec==None:
        prec=butools.checkPrecision

    return CheckRAPRepresentation(H[0],SumMatrixList(H[1:]),prec)
