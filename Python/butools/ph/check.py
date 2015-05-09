# -*- coding: utf-8 -*-
"""
Created on Sun Aug 24 15:31:04 2014

@author: gabor
"""

import butools
from butools.mc import CheckGenerator, CheckProbVector
import numpy as np
import scipy.linalg as la

def CheckPHRepresentation (alpha, A, prec=None):
    """
    Checks if the given vector and matrix define a valid phase-
    type representation.
    
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
        A is a transient generator, and they have the same size.
    """

    if prec==None:
        prec=butools.checkPrecision

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

def CheckMERepresentation (alpha, A, prec=None):
    """
    Checks if the given vector and matrix define a valid matrix-
    exponential representation.
    
    Parameters
    ----------
    alpha : matrix, shape (1,M)
        Initial vector of the matrix-exponential distribution 
        to check
    A : matrix, shape (M,M)
        Matrix parameter of the matrix-exponential distribution
        to check
    prec : double, optional
        Numerical precision. The default value is 1e-14.
    
    Returns
    -------
    r : bool
        True, if the matrix is a square matrix, the vector and 
        the matrix have the same size, the dominant eigenvalue
        is negative and real
    
    Notes
    -----
    This procedure does not check the positivity of the density!
    Call 'CheckMEPositiveDensity' if it is needed, but keep in
    mind that it can be time-consuming, while this procedure
    is fast.
    """

    if prec==None:
        prec=butools.checkPrecision

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
    
