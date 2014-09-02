# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 19:15:16 2013

@author: gabor
"""

import numpy as np
import numpy.matlib as ml
import numpy.linalg as la

def MStaircase (Y,Z, precision=1e-12):
    """
    Computes a smaller representation using the staircase 
    algorithm.
    
    Notes
    -----
    This function should not be called directly.
    It is used by 'MinimalRepFromME' and 'MinimalRepFromRAP'.
    
    References
    ----------
    .. [1]  P. Buchholz, M. Telek, "On minimal representation 
            of rational arrival processes." Madrid Conference
            on Qeueuing theory (MCQT), June 2010.
    """

    X = []
    for y in Y:
        X.append(ml.matrix(y))
    M = len(X)
    m = X[0].shape[0]
    U = ml.eye(m)

    #The sum of the ranks calculated in every loop
    ranksum = 0
    crit = True #The stopping criteria
    while crit==True:
        r = la.matrix_rank(Z, tol=precision)
        ranksum += r
        [Ui,S,T] = la.svd(Z)

        Transf = ml.eye(ranksum - r + Ui.shape[0])
        Transf[-Ui.shape[1]:,-Ui.shape[0]:] = Ui.T
        U = (Transf*U.T).T

        for i in range(M):
            TEMP = Ui.T*X[i]*Ui
            X[i] = TEMP[r:,r:]
            if i==0:
                Z = TEMP[r:,0:r]
            else:
                Z = np.hstack((Z,TEMP[r:,0:r]))

        if la.norm(Z) < precision or la.matrix_rank(Z, tol=precision) == m-ranksum:
            crit = False

    n = ranksum  
    if la.norm(Z) < precision:
        n = ranksum
        x = np.sum(U.T,1)
        x = x[0:n]
        
        #does x have a 0 value somewhere
        yes = False
        zeroloc = ml.zeros((n,1))
        nonzero = -1 # this will indicate a row of x for which x's value is non-zero
        for l in range(n):
           if abs(x[l]) < precision:
               yes = True
               zeroloc[l,0] = 1
           elif nonzero==-1:
               nonzero = l
        
        R = ml.eye(n)
        if yes:
            for l in range(n):
               if zeroloc[l,0]==1:
                   R[l,nonzero] = 1
        y = R*x
        Gamma = ml.diag(np.array(y).flatten())
        TEMP1 = ml.eye(m)
        TEMP1[0:n,0:n] = la.inv(Gamma)
        TEMP2 = ml.eye(m)
        TEMP2[0:n,0:n] = R
        
        B = la.inv (TEMP1*TEMP2*U.T)
    else:
        n=m
        B=ml.eye(m)
    return (B,n)
