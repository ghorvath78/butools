# -*- coding: utf-8 -*-
"""
Created on Sun Aug 24 15:31:04 2014

@author: gabor
"""

import butools
from butools.mc import CheckProbMatrix, CheckProbVector
import numpy as np
import scipy.linalg as la



def CheckDPHRepresentation (alpha, A, prec=None):
    """
    Checks if the given vector and matrix define a valid 
    discrete phase-type representation.
    
    Parameters
    ----------
    alpha : matrix, shape (1,M)
        Initial vector of the phase-type distribution to check
    A : matrix, shape (M,M)
        Transient generator of the phase-type distribution to
        check
    prec : double, optional
        Numerical precision. The default value is 1e-14.
    
    Returns
    -------
    r : bool
        True, if vector alpha is a probability vector and matrix
        A is substochastic, and they have the same size.
    """

    if prec==None:
        prec=butools.checkPrecision

    if len(alpha.shape)<2:
        if butools.verbose:
            print("CheckDPHRepresentation: Initial vector must be a matrix of size (1,N)!")
        return False
    
    if alpha.shape[1]!=A.shape[0]:
        if butools.verbose:
            print("CheckDPHRepresentation: The vector and the matrix have different sizes!")
        return False

    if not CheckProbMatrix(A,True,prec) or not CheckProbVector(alpha,True,prec):
        return False

    return True

def CheckMGRepresentation (alpha, A, prec=None):
    """
    Checks if the given vector and matrix define a valid matrix-
    geometric representation.
    
    Parameters
    ----------
    alpha : matrix, shape (1,M)
        Initial vector of the matrix-geometric distribution 
        to check
    A : matrix, shape (M,M)
        Matrix parameter of the matrix-geometric distribution
        to check
    prec : double, optional
        Numerical precision. The default value is 1e-14.
    
    Returns
    -------
    r : bool
        True, if the matrix is a square matrix, the vector and 
        the matrix have the same size, the dominant eigenvalue
        is positive, less than 1 and real. 
    
    Notes
    -----
    This procedure does not check the positivity of the density!
    The discrete counterpart of 'CheckMEPositiveDensity' does
    not exist yet (research is needed).
    """

    if prec==None:
        prec=butools.checkPrecision

    if len(alpha.shape)<2:
        if butools.verbose:
            print("CheckMGRepresentation: Initial vector must be a matrix of size (1,N)!")
        return False

    if A.shape[0]!=A.shape[1]:
        if butools.verbose:
            print("CheckMGRepresentation: The matrix is not a square matrix!")
        return False

    if alpha.shape[1]!=A.shape[0]:
        if butools.verbose:
            print("CheckMGRepresentation: The vector and the matrix have different sizes!")
        return False

    if np.sum(alpha)<-prec*alpha.size or np.sum(alpha)>1.0+prec*alpha.size:
        if butools.verbose:
            print("CheckMGRepresentation: The sum of the vector elements is less than zero or greater than one!")
        return False

    ev = la.eigvals(A)
    ix = np.argsort(-np.abs(np.real(ev)))
    maxev = ev[ix[0]]

    if not np.isreal(maxev):
        if butools.verbose:
            print("CheckMGRepresentation: The largest eigenvalue of the matrix is complex!")
        return False

    if maxev>1.0+prec:
        if butools.verbose:
            print("CheckMGRepresentation: The largest eigenvalue of the matrix is greater than 1!")
        return False       

    if np.sum(np.abs(ev)==abs(maxev)) > 1 and butools.verbose:
        print ("CheckMGRepresentation warning: There are more than one eigenvalue with the same absolute value as the largest eigenvalue!")

    return True
    
