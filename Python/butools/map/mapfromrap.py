# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 08:21:05 2013

@author: gabor
"""

import numpy as np
from numpy import linalg as la
import butools
from butools.reptrans import FindMarkovianRepresentation
from butools.map import CheckRAPRepresentation, CheckMRAPRepresentation

def MMAPFromMRAP (H, prec=1e-14):
    """
    Obtains a Markovian representation of a rational
    arrival process of the same size, if possible, using the
    procedure published in [1]_.
    
    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the MRAP to transform
    precision : double, optional
        A representation is considered to be a Markovian one
        if it is closer to it than this precision
    
    Returns
    -------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the MMAP (if found)
    
    References
    ----------
    .. [1] András Horváth, Gábor Horváth, Miklós Telek, "A 
           traffic based decomposition of two-class queueing 
           networks with priority service". COMPUTER NETWORKS 
           53:(8) pp. 1235-1248. (2009)
    """

    if butools.checkInput and not CheckMRAPRepresentation (H):
        raise Exception("MMAPFromMRAP: Input is not a valid MRAP representation!")    

    def transfun (oH, B):
        return [la.inv(B)*oHk*B for oHk in oH]
        
    def evalfun (oH, k=0):
        oH0 = oH[0] - np.diag(np.diag(oH[0]))
        if k%2 == 0:
            dist = np.min(oH0)
            for oHk in oH[1:]:
                dist = min(dist, np.min(oHk))
            return -dist
        else:
            dist = np.sum(oH0[oH0<0])
            for oHk in oH[1:]:
                dist += np.sum(oHk[oHk<0])
            return -dist

    return FindMarkovianRepresentation (H, transfun, evalfun, prec)

def MAPFromRAP (H0, H1, prec=1e-14):
    """
    Obtains a Markovian representation of a rational
    arrival process of the same size, if possible, using the
    procedure published in [1]_.
    
    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational arrival process
    precision : double, optional
        A representation is considered to be a Markovian one
        if it is closer to it than this precision
    
    Returns
    -------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process
    
    References
    ----------
    .. [1] G Horvath, M Telek, "A minimal representation of 
           Markov arrival processes and a moments matching 
           method," Performance Evaluation 64:(9-12) pp. 
           1153-1168. (2007)       
    """

    if butools.checkInput and not CheckRAPRepresentation (H0, H1):
        raise Exception("MAPFromRAP: Input is not a valid RAP representation!")    

    return MMAPFromMRAP([H0,H1], prec)
