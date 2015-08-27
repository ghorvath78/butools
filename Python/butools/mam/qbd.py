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
    """
    Returns the fundamental matrices corresponding to the
    given QBD. Matrices R, G and U can be returned
    depending on the "matrices" parameter.
    
    The implementation is based on [1]_, please cite it if
    you use this method.
    
    Parameters
    ----------
    B : matrix, shape (N,N)
        The matrix corresponding to backward transitions
    L : matrix, shape (N,N)
        The matrix corresponding to local transitions
    F : matrix, shape (N,N)
        The matrix corresponding to forward transitions
    matrices : string
        Specifies which matrices are required. 'R' means 
        that only matrix 'R' is needed. 'UG' means that
        matrices U and G are needed. Any combinations of
        'R', 'G' and 'U' are allowed, in any order.
    precision : double, optional
        The matrices are computed iteratively up to this
        precision. The default value is 1e-14
    maxNumIt : int, optional
        The maximal number of iterations. The default value
        is 50.
    method : {"CR", "LR", "NI", "FI", "IS"}, optional
        The method used to solve the matrix-quadratic
        equation (CR: cyclic reduction, LR: logarithmic
        reduction, NI: Newton iteration, FI: functional
        iteration, IS: invariant subspace method). The 
        default is "CR". "CR" is the only supported 
        method in the Mathematica and Python implementation.
    
    Returns
    -------
    M : list of matrices
        The list of calculated matrices in the order as
        requested in the 'matrices' parameter.
    
    Notes
    -----
    Discrete and continuous QBDs are both supported, the
    procedure auto-detects it based on the diagonal entries
    of matrix L.
    
    References
    ----------
    .. [1] Bini, D. A., Meini, B., Steffé, S., Van Houdt,
           B. (2006, October). Structured Markov chains 
           solver: software tools. In Proceeding from the
           2006 workshop on Tools for solving structured 
           Markov chains (p. 14). ACM.
    """

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
    """
    Returns the parameters of the matrix-geometrically 
    distributed stationary distribution of a QBD.
    
    Using vector pi0 and matrix R provided by this function
    the stationary solution can be obtained by
    
    .. math::
        \pi_k=\pi_0 R^k.    
    
    Parameters
    ----------
    B : matrix, shape (N,N)
        The matrix corresponding to backward transitions
    L : matrix, shape (N,N)
        The matrix corresponding to local transitions
    F : matrix, shape (N,N)
        The matrix corresponding to forward transitions
    L0 : matrix, shape (N,N)
        The matrix corresponding to local transitions at
        level zero
    precision : double, optional
        The fundamental matrix R is computed up to this
        precision. The default value is 1e-14
    
    Returns
    -------
    pi0 : matrix, shape (1,N)
        The stationary probability vector of level zero
    R : matrix, shape (N,N)
        The matrix parameter of the matrix geometrical
        distribution of the QBD 
    """
    
    m = L0.shape[0]
    I = ml.eye(m)
    
    R = QBDFundamentalMatrices (B, L, F, "R", prec)
    
    # Convert to discrete time problem, if needed
    if np.sum(np.diag(L0) < 0): # continues time
        lamb = float(np.max(-np.diag(L0)))
        B = ml.matrix(B) / lamb
        L0 = ml.matrix(L0) / lamb + I
    
    pi0 = DTMCSolve(L0+R*B)
    nr = np.sum(pi0*la.inv(I-R))
    pi0 /= nr
    
    return pi0, R

def QBDStationaryDistr (pi0, R, K):
    """
    Returns the stationary distribution of a QBD up to a
    given level K.
    
    Parameters
    ----------
    pi0 : matrix, shape (1,N)
        The stationary probability vector of level zero
    R : matrix, shape (N,N)
        The matrix parameter of the matrix geometrical
        distribution of the QBD 
    K : integer
        The stationary distribution is returned up to
        this level.
    
    Returns
    -------
    pi : array, length (K+1)*N
        The stationary probability vector up to level K
    """

    m = R.shape[0]    
    qld = ml.empty((1,(K+1)*m))
    qld[0,0:m] = pi0
    pix = pi0
    for k in range(1,K+1):
        pix = pix*R
        qld[0,k*m:(k+1)*m] = pix
    return qld
