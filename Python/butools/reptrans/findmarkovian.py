# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 08:24:03 2013

@author: gabor
"""
import copy
import numpy as np
import numpy.matlib as ml

def FindMarkovianRepresentation (rep, transfun, evalfun, precision=1e-7):
    """
    Obtains a Markovian representation from a non-Markovian 
    one while keeping the size the same, by applying a series 
    of elementary transformations.
    
    Parameters
    ----------
    rep : tuple of matrices
        The initial non-Markovian representation
        (initial vector and generator of a PH, matrices of a
        MAP, or a MMAP, etc.)
    transfun : callable
        A function that transforms the representation using 
        the given similarity transformation matrix
    evalfunc : callable
        A function that returns how far the representation is
        from the Markovian one
    precision : double
        A representation is considered to be a Markovian one
        if it is closer than the precision. The default value
        is 1e-7
        
    Returns
    -------
    mrep : tuple of matrices
        The Markovian representation, if found. If not found,
        the closest one is returned.
    
    Notes
    -----
    This function should not be called directly.
    It is used by 'PHFromME', 'MAPFromRAP', etc. functions.
    
    References
    ----------
    .. [1]  G Horváth, M Telek, "A minimal representation of 
            Markov arrival processes and a moments matching 
            method," Performance Evaluation 64:(9-12) 
            pp. 1153-1168. (2007)
    """

    def elementary (erep, b, k):
        bestdist = evalfun (erep, k)
        bestrep = erep
        repSize = erep[0].shape[1]
        for i in range(repSize):
            for j in range(repSize):
                if i!=j:
                    # create elementary transformation matrix with +b
                    B = ml.eye(repSize)
                    B[i,j] = b
                    B[i,i] = 1.0 - b
                    # apply similarity transform
                    newrep = transfun (erep, B)
                    newdist = evalfun (newrep, k)
                    # store result if better
                    if newdist < bestdist:
                        bestrep = newrep
                        bestdist = newdist
                    # create elementary transformation matrix with -b
                    B = ml.eye(repSize)
                    B[i,j] = -b
                    B[i,i] = 1.0 + b
                    # apply similarity transform
                    newrep = transfun (erep, B)
                    newdist = evalfun (newrep, k)
                    # store result if better
                    if newdist < bestdist:
                        bestrep = newrep
                        bestdist = newdist
        return (bestrep, bestdist)

    def minimize (orep, iters, b, k):
        lastdist = evalfun(orep, k)
        bestrep = orep
        for i in range(iters):
            orep, dist = elementary (orep, b, k)
            if dist >= lastdist:
                break
            else:
                lastdist = dist
                bestrep = orep
        return (bestrep, lastdist)

    if evalfun(rep) < precision:
        return rep
    
    nrep = copy.deepcopy(rep)
    M = nrep[0].shape[1]
    b = 0.5
    odist = np.inf
    while b>precision/2:
        for j in range(1,M*M):
            for k in range(4):
                nrep, dist = minimize(nrep, M*M, b, k)
                if dist < precision:
                    return nrep
            if odist <= dist:
                break
            odist = dist
        b = b / 2.0
    return nrep
