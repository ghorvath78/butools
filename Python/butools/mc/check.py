# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 16:57:48 2014

@author: gabor
"""

import numpy as np
import numpy.linalg as la
import butools

def CheckGenerator (Q, transient=False, prec=None):
    """
    Checks if the matrix is a valid generator matrix: the 
    matrix is a square matrix, the matrix has positive or 
    zero off-diagonal elements, the diagonal of the matrix 
    is negative, the rowsum of the matrix is 0.
    
    If the "transient" parameter is set to false, it checks 
    if the real part of the maximum absolute eigenvalue is 
    less than zero and the rowsum is equal or less than 0. 
    
    Parameters
    ----------
    Q : matrix, shape (M,M)
        The generator to check.
    transient : bool, optional
        If true, the procedure checks if Q is a transient 
        generator, otherwise it checks if it is a valid 
        generator. The default value is false.
    prec : double, optional
        Entries with absolute value less than prec are 
        considered to be zeros. The default value is 1e-14.
        
    Returns
    -------
    r : bool
        The result of the check.
    """

    if prec==None:
        prec=butools.checkPrecision

    if not isinstance(Q,np.ndarray):
        Q = np.array(Q)
        
    if Q.shape[0]!=Q.shape[1]:
        if butools.verbose:
            print ("CheckGenerator: Generator is not a square matrix!")
        return False

    if np.any(np.diag(Q)>=prec):
        if butools.verbose:
            print ("CheckGenerator: The diagonal of the generator is not negative (precision: {0})!".format(prec))
        return False

    N = Q.shape[0]
    odQ = Q<-prec
    for i in range(N):
        odQ[i,i] = 0

    if np.sum(np.any(odQ))>0:
        if butools.verbose:
            print ("CheckGenerator: The generator has negative off-diagonal element (precision: {0})!".format(prec))
        return False

    if transient:
        if np.max(np.sum(Q,1))>prec:
            if butools.verbose:
                print ("CheckGenerator: The rowsum of the transient generator is greater than 0 (precision: {0})!".format(prec))
            return False

        if np.max(np.real(la.eigvals(Q)))>=prec:
            if butools.verbose:
                print ("CheckGenerator: The transient generator has non-negative eigenvalue (precision: {0})!".format(prec))
            return False
    else:
        if np.any(np.abs(np.sum(Q,1))>prec):
            if butools.verbose:
                print ("CheckGenerator: The rowsum of the generator is not 0 (precision: {0})!".format(prec))
            return False
    return True

def CheckProbMatrix (P, transient=False, prec=None):
    """
    Checks if the matrix is a valid probability matrix: the 
    matrix is a square matrix, the matrix has positive or 
    zero off-diagonal elements, the rowsum of the matrix is 1.
    
    If "transient" is true, it checks if the matrix is a 
    valid transient probability matrix: the matrix is a square
    matrix, the matrix has positive or zero off-diagonal 
    elements, the rowsum of the matrix is less than or equal
    to 1, the maximum absolute eigenvalue is less than 1. 
    
    Parameters
    ----------
    P : matrix, shape (M,M)
        The matrix to check.
    transient : bool, optional
        If true, the procedure checks if P is a transient 
        probability matrix, otherwise it checks if it is
        a valid probability matrix. The default value is 
        false.
    prec : double, optional
        Entries with absolute value less than prec are 
        considered to be zeros. The default value is 1e-14.
        
    Returns
    -------
    r : bool
        The result of the check.
    """

    if prec==None:
        prec=butools.checkPrecision

    if not isinstance(P,np.ndarray):
        P = np.array(P)

    if P.shape[0]!= P.shape[1]:
        if butools.verbose:
            print ("CheckProbMatrix: the matrix is not a square matrix!")
        return False

    if np.min(np.min(P))<-prec:
        if butools.verbose:
            print ("CheckProbMatrix: the matrix has negative element (precision: {0})!".format(prec))
        return False

    if transient:
        if np.any(np.sum(P,1)-1.0>prec*P.shape[1]):
            if butools.verbose:
                print ("CheckProbMatrix: The rowsum of the matrix (transient) is not less or equal than 1 (precision: {0})!", prec)
            return False

        if np.max(np.real(la.eigvals(P)))>=1.0-prec:
            if butools.verbose:
                print ("CheckProbMatrix: The real part of the largest eigenvalue of the transient matrix is not less than 1 (precision: {0})!".format(prec))
            return False
    else:
        if np.any(np.abs(np.sum(P,1)-1.0)>prec*P.shape[1]):
            if butools.verbose:
                print ("CheckProbMatrix: The rowsum of the matrix is not 1 (precision: {0})!".format(prec))
            return False
    return True

def CheckProbVector (pi, sub=False, prec=None):
    """
    Checks if the vector is a valid probability vector: the 
    vector has only non-negative elements, the sum of the 
    vector elements is 1.
    
    If parameter "sub" is set to true, it checks if the 
    vector is a valid substochastic vector: the vector has 
    only non-negative elements, the sum of the elements are
    less than 1.
    
    Parameters
    ----------
    pi : vector, shape (1, M) or (M, 1)
        The matrix to check.
    sub : bool, optional
        If false, the procedure checks for stochastic, if 
        true, it checks for sub-stochastic property. The 
        default value is false.
    prec : double, optional
        Numerical precision. Entries with absolute value 
        less than prec are considered to be zeros. The 
        default value is 1e-14.
        
    Returns
    -------
    r : bool
        The result of the check.
    """

    if prec==None:
        prec=butools.checkPrecision

    if np.min(pi)<-prec:
        if butools.verbose:
            print ("CheckProbVector: The vector has negative element (precision: {0})!".format(prec))
        return False

    if sub:
        if np.sum(pi)>1.0+prec*pi.size:
            if butools.verbose:
                print ("CheckProbVector: The sum of the substochastic vector is not less than 1 (precision: {0})!".format(prec))
            return False
    else:
        if np.abs(np.sum(pi)-1.0)>prec*pi.size:
            if butools.verbose:
                print ("CheckProbVector: The sum of the vector is not 1 (precision: {0})!".format(prec))
            return False
    
    return True
