# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 15:19:44 2013

@author: gabor
"""

import numpy as np
import scipy.linalg as la
import numpy.linalg as nla
from numpy import matlib as ml
from butools.moments import ReducedMomsFromMoms
from butools.reptrans import MStaircase
import butools
from butools.ph import MomentsFromME, MEFromMoments, CheckMERepresentation

def MEOrderFromMoments (moms, prec=1e-12):
    """
    Returns the order of ME distribution that can realize
    the given moments.
    
    Parameters
    ----------
    moms : list of doubles
        The list of moments
    prec : double, optional
        Precision used to detect if the determinant of the
        Hankel matrix is zero. The default value is 1e-12.
    
    Returns
    -------
    order : int
        The order of ME distribution that can realize the 
        given moments
    
    References
    ----------
    .. [1]  L. Bodrog, A. Horvath, M. Telek, "Moment 
            characterization of matrix exponential and Markovian
            arrival processes," Annals of Operations Research, 
            vol. 160, pp. 51-68, 2008.
    """

    sizem=int((len(moms)+1)/2)
    rmoms=[1] + ReducedMomsFromMoms(moms);
    for k in range(1,sizem+1):
        hankel = np.zeros((k,k))
        for i in range(k):
           for j in range(k):
               hankel[i,j] = rmoms[i+j]
        if abs(la.det(hankel)) < prec:
           return k-1
    return sizem
        
def MEOrder (alpha, A, kind="moment", prec=1e-10):
    """
    Returns the order of the ME distribution (which is not 
    necessarily equal to the size of the representation).
    
    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-exponential 
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential 
        distribution.
    kind : {'obs', 'cont', 'obscont', 'moment'}, optional
        Determines which order is computed. Possibilities: 
        'obs': observability, 
        'cont': controllability,
        'obscont': the minimum of observability and 
        controllability order,
        'moment': moment order (which is the default).
    prec : double, optional
        Precision used to detect if the determinant of the 
        Hankel matrix is zero (in case of kind="moment" only),
        or the tolerance for the rank calculation. The
        default value is 1e-10.
    
    Returns
    -------
    order : int
        The order of ME distribution
    
    References
    ----------
    .. [1]  P. Buchholz, M. Telek, "On minimal representation
            of rational arrival processes." Madrid Conference on
            Qeueuing theory (MCQT), June 2010.
    """

    if butools.checkInput and not CheckMERepresentation (alpha, A):
        raise Exception("MEOrder: Input is not a valid ME representation!")

    N = alpha.shape[1]
    if kind=="cont":
        re = np.zeros ((N,N))
        for n in range(N):
            re[n,:] = np.sum(A.T**n, 0)
        return nla.matrix_rank (re, prec)
    elif kind=="obs":
        re = np.zeros ((N,N))
        for n in range(N):
            re[n,:] = alpha*A**n
        return nla.matrix_rank (re, prec)
    elif kind=="obscont":
        re = np.zeros ((N,N))
        for n in range(N):
            re[n,:] = alpha*A**n
        obsOrder = nla.matrix_rank (re, prec)
        re = np.zeros ((N,N))
        for n in range(N):
            re[n,:] = np.sum(A.T**n, 0)
        contOrder = nla.matrix_rank (re, prec)
        return min(obsOrder,contOrder)        
    elif kind=="moment":
        return MEOrderFromMoments (MomentsFromME(alpha, A), prec)
    else:
        raise Exception("Invalid 'kind' parameter!")

def MinimalRepFromME (alpha, A, how="moment", precision=1e-12):
    """
    Returns the minimal representation of the given ME 
    distribution.
    
    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-exponential 
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential 
        distribution.
    how : {"obs", "cont", "obscont", "moment"}, optional        
        Determines how the representation is minimized. 
        Possibilities:
        'obs': observability, 
        'cont': controllability,
        'obscont': the minimum of observability and 
        controllability order,
        'moment': moment order (which is the default).
    precision : double, optional
       Precision used by the Staircase algorithm. The default
       value is 1e-12.
    
    Returns
    -------
    beta : vector, shape (1,N)
        The initial vector of the minimal representation
    B : matrix, shape (N,N)
        The matrix parameter of the minimal representation
    
    References
    ----------
    .. [1]  P. Buchholz, M. Telek, "On minimal representation
            of rational arrival processes." Madrid Conference on
            Qeueuing theory (MCQT), June 2010.
    """

    if butools.checkInput and not CheckMERepresentation (alpha, A):
        raise Exception("MinimalRepFromME: Input is not a valid ME representation!")

    if how=="cont":
        H0 = A
        H1 = np.sum(-A,1) * alpha
        B, n = MStaircase ([H0, H1], ml.ones((A.shape[0],1)), precision)
        return ((alpha*B)[0,0:n], (la.inv(B)*A*B)[0:n,0:n])
    elif how=="obs":
        H0 = A
        H1 = np.sum(-A,1) * alpha
        G = [H0.T,H1.T]
        B, n = MStaircase (G, alpha.T, precision)
        return ((alpha*B)[:,0:n], (la.inv(B)*A*B)[0:n,0:n])
    elif how=="obscont":
        alphav, Av = MinimalRepFromME (alpha, A, "cont", precision)
        return MinimalRepFromME (alphav, Av, "obs", precision)      
    elif how=="moment":
        N = MEOrder (alpha, A, "moment", precision)
        moms = MomentsFromME (alpha, A, 2*N-1)
        return MEFromMoments (moms)
        