# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 18:29:52 2013

@author: gabor
"""
import numpy as np
import numpy.matlib as ml
import numpy.linalg as la
import butools
from butools.mc import CheckGenerator, CheckProbMatrix


def CRPSolve (Q):
    """
    Computes the stationary solution of a continuous time 
    rational process (CRP).
    
    Parameters
    ----------
    Q : matrix, shape (M,M)
        The generator matrix of the rational process
        
    Returns
    -------
    pi : row vector, shape (1,M)
        The vector that satisfies 
        `\pi\, Q = 0, \sum_i \pi_i=1`
    
    Notes
    -----
    Continuous time rational processes are like continuous 
    time Markov chains, but the generator does not have to 
    pass the :func:`CheckGenerator` test (but the rowsums 
    still have to be zeros).
    """
    
    if butools.checkInput and np.any(np.sum(Q,1)>butools.checkPrecision):
        raise Exception("CRPSolve: The matrix has a rowsum which isn't zero!")
    
    M = np.array(Q)
    M[:,0] = np.ones(M.shape[0])
    m = np.zeros(M.shape[0])
    m[0] = 1.0
    return ml.matrix(la.solve (M.T, m))
#    return la.solve (M.T, m)
    
def CTMCSolve (Q):
    """
    Computes the stationary solution of a continuous time 
    Markov chain.
    
    Parameters
    ----------
    Q : matrix, shape (M,M)
        The generator matrix of the Markov chain
        
    Returns
    -------
    pi : row vector, shape (1,M)
        The vector that satisfies `\pi\, Q = 0, \sum_i \pi_i=1`
    
    Notes
    -----
    The procedure raises an exception if :code:`checkInput` 
    is set to :code:`true` and :func:`CheckGenerator` (Q) fails.
    """

    if butools.checkInput and not CheckGenerator(Q, False):
        raise Exception("CTMCSolve: The given matrix is not a valid generator. If you are sure you want this use CRPSolve instead of CTMCSolve.")

    return CRPSolve(Q)
    
def DRPSolve (P):
    """
    Computes the stationary solution of a discrete time 
    Markov chain.
    
    Parameters
    ----------
    P : matrix, shape (M,M)
        The matrix parameter of the rational process
        
    Returns
    -------
    pi : row vector, shape (1,M)
        The vector that satisfies 
        `\pi\, P = \pi, \sum_i \pi_i=1`
    
    Notes
    -----
    Discrete time rational processes are like discrete time 
    Markov chains, but the P matrix does not have to pass 
    the :func:`CheckProbMatrix` test (but the rowsums still 
    have to be ones).
    """

    if butools.checkInput and np.any(np.sum(P,1)-1.0>butools.checkPrecision):
        raise Exception("DRPSolve: The matrix has a rowsum which isn't 1!")

    if not isinstance(P,np.ndarray):
        P = np.array(P)

    return CRPSolve(P-ml.eye(P.shape[0]))
    
def DTMCSolve (P):
    """
    Computes the stationary solution of a discrete time 
    Markov chain.
    
    Parameters
    ----------
    P : matrix, shape (M,M)
        The transition probability matrix of the Markov 
        chain
        
    Returns
    -------
    pi : row vector, shape (1,M)
        The vector that satisfies `\pi\, P = \pi, \sum_i \pi_i=1`
    
    Notes
    -----
    The procedure raises an exception if :code:`butools.checkInput` 
    is set to :code:`true` and :func:`CheckProbMatrix` (P) fails.
    """

    if butools.checkInput and not CheckProbMatrix(P, False):
        raise Exception("CTMCSolve: The given matrix is not a valid generator. If you are sure you want this use CRPSolve instead of CTMCSolve.")

    return DRPSolve(P)

