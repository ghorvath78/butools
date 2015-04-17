# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 17:44:47 2013

@author: gabor
"""

import numpy as np
import numpy.matlib as ml
import numpy.linalg as la
import butools
from butools.mc import DTMCSolve

def QBDFundamentalMatrices (B, L, F, matrices="G", precision=1e-14, maxNumIt=50, method="CR", shift=True):

    m = L.shape[0]
    I = ml.eye(m)
    
    if method!="CR" and butools.verbose:
        print("Warning: Currently only the 'CR' method is available in the Python implementation!")        

    # Convert to discrete time problem, if needed
    continuous = False
    if np.sum(np.diag(L) < 0): # continues time
        continuous = True
        lamb = float(np.max(-np.diag(L)))
        B = ml.matrix(B) / lamb
        L = ml.matrix(L) / lamb + I
        F = ml.matrix(F) / lamb
    
    # Shift technique
    if shift:
        theta = DTMCSolve (B+L+F)
        drift = theta * np.sum(B,1) - theta * np.sum(F,1)
        Lold = L
        if drift < 0: # MC is transient -> use the dual MC
            Fold = F
            F = F - ml.ones((m,1))*(theta*F)
            L = L + ml.ones((m,1))*(theta*B)
        else:
            uT = ml.ones((1,m)) / m
            Bold = B
            B = B - np.sum(B,1)*uT
            L = L + np.sum(F,1)*uT

    # Start of Logaritmic Reduction (Basic)
    BF = la.inv(I - L)
    BB = BF * F
    BF = BF * B
    G = BF
    PI = BB
    check = 1.0
    numit = 0
    while check > precision and numit < maxNumIt:
        Lstar = BF*BB + BB*BF
        Bstar = BB*BB
        Fstar = BF*BF
        BB = la.inv(I-Lstar)
        BF = BB * Fstar
        BB *= Bstar
        G += PI * BF
        PI *= BB
        check = min(la.norm(BB,np.inf),la.norm(BF,np.inf))
        numit += 1


    if numit == maxNumIt and butools.verbose==True:
        print("Maximum Number of Iterations reached")

    # Shift Technique
    if shift==True:
        L = Lold
        if drift < 0: # transient
            F = Fold  # restore original A2
        else:  # pos recurrent
            G = G + ml.ones((m,1))*uT
            B = Bold  # restore original A0

    if butools.verbose==True:
        res_norm = la.norm (G-B-(L+F*G)*G, np.inf)
        print("Final Residual Error for G: ", res_norm)

    ret = []
    for M in matrices:
        if M=="G":
            ret.append(G)
        elif M=="R":
            R = F*la.inv(I-(L+F*G))
            if butools.verbose==True:
                res_norm = la.norm (R-F-R*(L+R*B), np.inf)
                print("Final Residual Error for R: ", res_norm)
            ret.append(R)
        elif M=="U":
            U = L+F*G
            if butools.verbose==True:
                res_norm = la.norm (U-L-F*la.inv(I-U)*B, np.inf)
                print("Final Residual Error for U: ", res_norm)
            if continuous==True:
                U = lamb*(U-I)
            ret.append(U)
    if len(ret)==1:
        return ret[0]
    else:
        return ret                

def QBDSolve (B, L, F, L0, prec=1e-14):
    
    m = L0.shape[0]
    I = ml.eye(m)
    
    R = QBDFundamentalMatrices (B, L, F, "R", prec)
    
    # Convert to discrete time problem, if needed
    if np.sum(np.diag(L0) < 0): # continues time
        lamb = float(np.max(-np.diag(L0)))
        B = ml.matrix(B) / lamb
        L0 = ml.matrix(L0) / lamb + I
    
    pi0 = DTMCSolve(L0+R*B, prec)
    nr = np.sum(pi0*la.inv(I-R))
    pi0 /= nr
    
    return pi0, R

def QBDStationaryDistr (pi0, R, K):

    m = R.shape[0]    
    qld = ml.empty((1,(K+1)*m))
    qld[0,0:m] = pi0
    pix = pi0
    for k in range(1,K+1):
        pix = pix*R
        qld[0,k*m:(k+1)*m] = pix
    return qld
