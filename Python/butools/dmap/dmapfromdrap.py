# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 08:21:05 2013

@author: gabor
"""

import numpy as np
from numpy import linalg as la
import numpy.matlib as ml
import butools
from butools.reptrans import FindMarkovianRepresentation
from butools.dmap import CheckDRAPRepresentation, CheckDMRAPRepresentation

def DMMAPFromDMRAP (H, prec=1e-14):
    """
    Obtains a Markovian representation of a discrete rational
    arrival process of the same size, if possible, using the
    procedure published in [1]_.
    
    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the DMRAP to transform
    precision : double, optional
        A representation is considered to be a Markovian one
        if it is closer to it than this precision
    
    Returns
    -------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the DMMAP (if found)
    
    References
    ----------
    .. [1] András Horváth, Gábor Horváth, Miklós Telek, "A 
           traffic based decomposition of two-class queueing 
           networks with priority service". COMPUTER NETWORKS 
           53:(8) pp. 1235-1248. (2009)
    """

    if butools.checkInput and not CheckDMRAPRepresentation (H):
        raise Exception("DMMAPFromDMRAP: Input is not a valid DMRAP representation!")    

    def transfun (oH, B):
        return [la.inv(B)*oHk*B for oHk in oH]
        
    def evalfun (oH, k=0):
        Ones = ml.ones(oH[0].shape)
        if k%2 == 0:
            dist = np.min(oH[0])
            for oHk in oH:
                dist = min(dist, np.min(oHk), np.min(Ones-oHk))
            return -dist
        else:
            dist = np.sum(oH[0][oH[0]<0])
            for oHk in oH:
                dist += min(np.sum(oHk[oHk<0]), np.sum(oHk[Ones-oHk<0]))
            return -dist

    return FindMarkovianRepresentation (H, transfun, evalfun, prec)

def DMAPFromDRAP (H0, H1, prec=1e-14):
    """
    Obtains a Markovian representation of a discrete rational
    arrival process of the same size, if possible, using the
    procedure published in [1]_.
    
    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the discrete rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the discrete rational arrival process
    precision : double, optional
        A representation is considered to be a Markovian one
        if it is closer to it than this precision
    
    Returns
    -------
    D0 : matrix, shape (M,M)
        The D0 matrix of the discrete Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the discrete Markovian arrival process
    
    References
    ----------
    .. [1] G Horvath, M Telek, "A minimal representation of 
           Markov arrival processes and a moments matching 
           method," Performance Evaluation 64:(9-12) pp. 
           1153-1168. (2007)       
    """

    if butools.checkInput and not CheckDRAPRepresentation (H0, H1):
        raise Exception("DMAPFromDRAP: Input is not a valid DRAP representation!")    

    return DMMAPFromDMRAP([H0,H1], prec)
