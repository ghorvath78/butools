# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 15:19:44 2013

@author: gabor
"""

import numpy as np
import scipy.linalg as la
from numpy import matlib as ml
from butools.moments import ReducedMomsFromMoms
from butools.reptrans import MStaircase
from butools.ph import MomentsFromME, MEFromMoments

def MEOrderFromMoments (moms, prec=1e-12):

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

    N = alpha.shape[1]
    if kind=="cont":
        re = np.zeros ((N,N))
        for n in range(N):
            re[n,:] = np.sum(A.T**n, 0)
        return la.matrix_rank (re, prec)
    elif kind=="obs":
        re = np.zeros ((N,N))
        for n in range(N):
            re[n,:] = alpha*A**n
        return la.matrix_rank (re, prec)
    elif kind=="obscont":
        re = np.zeros ((N,N))
        for n in range(N):
            re[n,:] = alpha*A**n
        obsOrder = la.matrix_rank (re, prec)
        re = np.zeros ((N,N))
        for n in range(N):
            re[n,:] = np.sum(A.T**n, 0)
        contOrder = la.matrix_rank (re, prec)
        return min(obsOrder,contOrder)        
    elif kind=="moment":
        return MEOrderFromMoments (MomentsFromME(alpha, A), prec)
    else:
        raise Exception("Invalid 'kind' parameter!")

def MinimalRepFromME (alpha, A, how="moment", precision=1e-12):

    if how=="cont":
        H0 = A
        H1 = np.sum(-A,1) * alpha
        B, n = MStaircase ([H0, H1], ml.ones((A.shape[0],1)), precision)
        return ((alpha*B)[0:n], (la.inv(B)*A*B)[0:n,0:n])
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
        moms = MomentsFromME (alpha, A, 2*N-1, precision)
        return MEFromMoments (moms)
        