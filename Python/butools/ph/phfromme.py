# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 08:21:05 2013

@author: gabor
"""

import numpy as np
from numpy import linalg as la
import butools
from butools.ph import CheckMERepresentation
from butools.reptrans import FindMarkovianRepresentation

def PHFromME (alpha, A, precision=1e-14):
    """
    Obtains a Markovian representation of a matrix 
    exponential distribution of the same size, if possible.
    
    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-exponential
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential 
        distribution.
    precision : double, optional
        A representation is considered to be a Markovian one
        if it is closer than the precision. The default value
        is 1e-14.
    
    Returns
    -------
    beta : vector, shape (1,M)
        The initial probability vector of the Markovian 
        monocyclic representation
    B : matrix, shape (M,M)
        Transient generator matrix of the Markovian 
        monocyclic representation
    
    References
    ----------
    .. [1] G Horváth, M Telek, "A minimal representation of 
           Markov arrival processes and a moments matching 
           method," Performance Evaluation 64:(9-12) pp. 
           1153-1168. (2007)
    """

    def transfun (orep, B):
        ao, Ao = orep
        return (ao*B, la.inv(B)*Ao*B)
        
    def evalfun (orep, k=0):
        ao, Ao = orep
        av= np.sum(-Ao,1)
        Ad = Ao-np.diag(np.diag(Ao))
        if k%2 == 0:
            return -min(np.min(ao), np.min(av), np.min(Ad))
        else:
            return -np.sum(ao[ao<0]) - np.sum(av[av<0]) - np.sum(Ad[Ad<0])

    if butools.checkInput and not CheckMERepresentation (alpha, A):
        raise Exception("PHFromME: Input is not a valid ME representation!")

    return FindMarkovianRepresentation ((alpha,A), transfun, evalfun, precision)
