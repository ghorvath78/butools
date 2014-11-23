# -*- coding: utf-8 -*-
"""
Created on Sun Aug 31 17:33:24 2014

@author: gabor
"""

import butools
import numpy as np
import numpy.matlib as ml
from butools.map import CheckMAPRepresentation, CheckMMAPRepresentation, MarginalDistributionFromMMAP
from numpy.random import rand
from butools.mc import CTMCSolve
from butools.utils import SumMatrixList


def SamplesFromMMAP (D, k, initial=None, prec=1e-14):
    """
    Generates random samples from a marked Markovian 
    arrival process.
    
    Parameters
    ----------
    D : list of matrices of shape(M,M), length(N)
        The D0...DN matrices of the MMAP
    K : integer
        The number of samples to generate.
    prec : double, optional
        Numerical precision to check if the input MMAP is
        valid. The default value is 1e-14.
    
    Returns
    -------
    x : matrix, shape(K,2)
        The random samples. Each row consists of two 
        columns: the inter-arrival time and the type of the
        arrival.        
    """

    if butools.checkInput and not CheckMMAPRepresentation (D, prec):
        raise Exception("SamplesFromMMAP: Input is not a valid MMAP representation!")    

    N = D[0].shape[0]
    
    if initial==None:
        # draw initial state according to the stationary distribution
        stst = CTMCSolve(SumMatrixList(D)).A.flatten()
        cummInitial = np.cumsum(stst)
        r = rand()
        state = 0
        while cummInitial[state]<=r:
            state+=1
    else:
        state = initial

    # auxilary variables
    sojourn = -1.0/np.diag(D[0])
    nextpr = ml.matrix(np.diag(sojourn))*D[0]
    nextpr = nextpr - ml.matrix(np.diag(np.diag(nextpr)))
    for i in range(1,len(D)):
        nextpr = np.hstack((nextpr, np.diag(sojourn)*D[i]))
    nextpr = np.cumsum(nextpr,1)
    
    if len(D)>2:
        x = np.empty((k,2))
    else:
        x = np.empty(k)

    for n in range(k):
        time = 0

        # play state transitions
        while state<N :
            time -= np.log(rand()) * sojourn[state]
            r = rand()
            nstate = 0
            while nextpr[state,nstate]<=r:
                nstate += 1
            state = nstate
        if len(D)>2:
            x[n,0] = time
            x[n,1] = state//N
        else:
            x[n] = time
        state = state % N
    
    return x

def SamplesFromMAP (D0, D1, k, initial=None, prec=1e-14):
    """
    Generates random samples from a Markovian arrival 
    process.
    
    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process
    K : integer
        The number of samples to generate.
    prec : double, optional
        Numerical precision to check if the input Markovian
        arrival process is valid. The default value is 
        1e-14.
    
    Returns
    -------
    x : vector, length(K)
        The vector of random samples (inter-arrival times).
    """

    if butools.checkInput and not CheckMAPRepresentation (D0, D1, prec):
        raise Exception("SamplesFromMAP: Input is not a valid MAP representation!")    

    return SamplesFromMMAP((D0,D1),k,initial,prec);
