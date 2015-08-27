# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 16:30:53 2013

@author: gabor
"""
import numpy as np
import scipy.linalg as la
import numpy.matlib as ml
import math
import butools
from butools.mc import CTMCSolve
from butools.utils import Diag, Linsolve

def FluidFundamentalMatrices (Fpp, Fpm, Fmp, Fmm, matrices, precision=1e-14, maxNumIt=50, method="CR"):
    """
    Returns the fundamental matrices corresponding to the
    given canonical Markov fluid model. Matrices Psi, K and
    U are returned depending on the "matrices" parameter.
    The canonical Markov fluid model is defined by the 
    matrix blocks of the generator of the background Markov
    chain partitioned according to the sign of the 
    associated fluid rates (i.e., there are "+" and "-" states).
    
    Parameters
    ----------
    Fpp : matrix, shape (Np,Np)
        The matrix of transition rates between states 
        having positive fluid rates
    Fpm : matrix, shape (Np,Nm)
        The matrix of transition rates where the source
        state has a positive, the destination has a 
        negative fluid rate associated.
    Fpm : matrix, shape (Nm,Np)
        The matrix of transition rates where the source
        state has a negative, the destination has a 
        positive fluid rate associated.
    Fpp : matrix, shape (Nm,Nm)
        The matrix of transition rates between states 
        having negative fluid rates
    matrices : string
        Specifies which matrices are required. 'P' means 
        that only matrix Psi is needed. 'UK' means that
        matrices U and K are needed. Any combinations of
        'P', 'K' and 'U' are allowed, in any order.
    precision : double, optional
        The matrices are computed iteratively up to this
        precision. The default value is 1e-14
    maxNumIt : int, optional
        The maximal number of iterations. The default value
        is 50.
    method : {"CR", "ADDA", "SDA"}, optional
        The method used to solve the algebraic Riccati
        equation (CR: cyclic reduction, ADDA: alternating-
        directional doubling algorithm, SDA: structured
        doubling algorithm). The default is "CR".
    
    Returns
    -------
    M : list of matrices
        The list of calculated matrices in the order as
        requested in the 'matrices' parameter.
    
    Notes
    -----
    Thanks to Benny Van Houdt for the implementation of the
    Riccati solvers.
    """

    if method=="CR":
        F = np.bmat ([[Fpp,Fpm],[Fmp,Fmm]])
        sp = Fpp.shape[0]
        sm = Fmm.shape[0]
        I = np.eye(sm+sp)
        Pf = I + F / np.max(-Diag(F))
        A1 = np.bmat([[Pf[:sp,:sp]/2.0, np.zeros((sp,sm))],[Pf[sp:,:sp], np.zeros((sm,sm))]])
        # Optimized cyclic reduction
        A1hat=ml.matrix(A1)
        R1 = np.eye(sp) / 2.0
        N = np.vstack((Pf[:sp,sp:]/2.0, Pf[sp:sp+sm,sp:sp+sm]))
        numit = 0
        while min(la.norm(R1,np.inf),la.norm(N,np.inf))>precision and numit<maxNumIt:
            Xmat = la.inv(I-A1)
            B = Xmat.dot(N)
            S1 = Xmat[sp:,:sp].dot(R1)
            U = R1.dot(B[:sp,:])
            A1 += np.hstack((N.dot(S1), np.vstack((U, np.zeros((sm,sm))))))
            A1hat[:sp,sp:] += U
            N = N.dot(B[sp:,:])
            R1 = R1.dot(Xmat[:sp,:sp]).dot(R1)
            numit += 1
        Ghat = la.inv(I-A1hat).dot(np.vstack((Pf[:sp,sp:]/2.0, Pf[sp:,sp:])))
        Psi = ml.matrix(Ghat[:sp,:])
    elif method=="ADDA" or method=="SDA":
        # via ADDA algorithm (Wang, Wang, Li 2011)
        A = ml.matrix(-Fpp)
        B = Fpm
        C = Fmp
        D = ml.matrix(-Fmm)
        gamma1 = np.max(np.diag(A))
        gamma2 = np.max(np.diag(D))
        if method=="SDA":
            gamma1 = max(gamma1,gamma2)
            gamma2 = gamma1
        sA = A.shape[0]
        sD = D.shape[0]
        IA = ml.eye(sA)
        ID = ml.eye(sD)
        A += gamma2*IA
        D += gamma1*ID
        Dginv = D.I
        Vginv = (D-C*A.I*B).I
        Wginv = (A-B*Dginv*C).I
        Eg = ID - (gamma1+gamma2)*Vginv
        Fg = IA - (gamma1+gamma2)*Wginv
        Gg = (gamma1+gamma2)*Dginv*C*Wginv
        Hg = (gamma1+gamma2)*Wginv*B*Dginv
        diff = 1.0
        numit = 0
        while diff>precision and numit<maxNumIt:
            Vginv = Eg*(ID-Gg*Hg).I
            Wginv = Fg*(IA-Hg*Gg).I
            Gg += Vginv*Gg*Fg
            Hg += Wginv*Hg*Eg
            Eg = Vginv*Eg
            Fg = Wginv*Fg
            neg = la.norm(Eg,1)
            nfg = la.norm(Fg,1)
            if method=="ADDA":
                eta = math.sqrt(nfg/neg)
                Eg *= eta
                Fg /= eta
                diff = neg*nfg
            else:
                diff = min(neg,nfg)
            numit += 1
        Psi = ml.matrix(Hg)
        
    if numit == maxNumIt and butools.verbose==True:
        print("Maximum Number of Iterations reached")
        
    if butools.verbose==True:
        res_norm = la.norm (Fpm+Fpp*Psi+Psi*Fmm+Psi*Fmp*Psi, np.inf)
        print("Final Residual Error for G: ", res_norm)

    ret = []
    for M in matrices:
        if M=="P":
            ret.append(Psi)
        elif M=="K":
            ret.append(Fpp+Psi*Fmp)
        elif M=="U":
            ret.append(Fmm+Fmp*Psi)
    if len(ret)==1:
        return ret[0]
    else:
        return ret                   
    
def FluidSolve (Fpp, Fpm, Fmp, Fmm, prec=1e-14):
    """
    Returns the parameters of the matrix-exponentially 
    distributed stationary distribution of a canonical 
    Markovian fluid model.
    
    The canonical Markov fluid model is defined by the 
    matrix blocks of the generator of the background Markov
    chain partitioned according to the sign of the 
    associated fluid rates (i.e., there are "+" and "-" states).   
    
    Using the returned 4 parameters the stationary
    solution can be obtained as follows.
    
    The probability that the fluid level is zero while 
    being in different states of the background process
    is given by vector mass0.
    
    The density that the fluid level is x while being in
    different states of the background process is
    
    .. math::
        \pi(x)=ini\cdot e^{K x}\cdot clo.    
    
    Parameters
    ----------
    Fpp : matrix, shape (Np,Np)
        The matrix of transition rates between states 
        having positive fluid rates
    Fpm : matrix, shape (Np,Nm)
        The matrix of transition rates where the source
        state has a positive, the destination has a 
        negative fluid rate associated.
    Fpm : matrix, shape (Nm,Np)
        The matrix of transition rates where the source
        state has a negative, the destination has a 
        positive fluid rate associated.
    Fpp : matrix, shape (Nm,Nm)
        The matrix of transition rates between states 
        having negative fluid rates
    precision : double, optional
        Numerical precision for computing the fundamental
        matrix. The default value is 1e-14
    
    Returns
    -------
    mass0 : matrix, shape (1,Np+Nm)
        The stationary probability vector of zero level
    ini : matrix, shape (1,Np)
        The initial vector of the stationary density
    K : matrix, shape (Np,Np)
        The matrix parameter of the stationary density
    clo : matrix, shape (Np,Np+Nm)
        The closing matrix of the stationary density
    """
    
    Psi, K, U = FluidFundamentalMatrices (Fpp, Fpm, Fmp, Fmm, "PKU", prec)
    
    mass0 = CTMCSolve(U)
    nr = np.sum(mass0) + 2.0*np.sum(mass0*Fmp*-K.I)
    
    return ml.hstack((ml.zeros((1,Fpp.shape[0])),mass0/nr)), mass0*Fmp/nr, K, ml.hstack((ml.eye(Fpp.shape[0]), Psi))

def GeneralFluidSolve (Q, R, Q0=[], prec=1e-14):
    """
    Returns the parameters of the matrix-exponentially 
    distributed stationary distribution of a general 
    Markovian fluid model, where the fluid rates associated
    with the states of the background process can be
    arbitrary (zero is allowed as well).
    
    Using the returned 4 parameters the stationary
    solution can be obtained as follows.
    
    The probability that the fluid level is zero while 
    being in different states of the background process
    is given by vector mass0.
    
    The density that the fluid level is x while being in
    different states of the background process is
    
    .. math::
        \pi(x)=ini\cdot e^{K x}\cdot clo.    
    
    Parameters
    ----------
    Q : matrix, shape (N,N)
        The generator of the background Markov chain
    R : diagonal matrix, shape (N,N)
        The diagonal matrix of the fluid rates associated
        with the different states of the background process
    Q0 : matrix, shape (N,N), optional
        The generator of the background Markov chain at 
        level 0. If not provided, or empty, then Q0=Q is 
        assumed. The default value is empty.
    precision : double, optional
        Numerical precision for computing the fundamental
        matrix. The default value is 1e-14
    
    Returns
    -------
    mass0 : matrix, shape (1,Np+Nm)
        The stationary probability vector of zero level
    ini : matrix, shape (1,Np)
        The initial vector of the stationary density
    K : matrix, shape (Np,Np)
        The matrix parameter of the stationary density
    clo : matrix, shape (Np,Np+Nm)
        The closing matrix of the stationary density
    """
    
    N = Q.shape[0]
    # partition the state space according to zero, positive and negative fluid rates
    ix = np.arange(N)
    ixz = ix[np.abs(np.diag(R))<=prec]
    ixp = ix[np.diag(R)>prec]
    ixn = ix[np.diag(R)<-prec]
    Nz = len(ixz)
    Np = len(ixp)
    Nn = len(ixn)
    # permutation matrix that converts between the original and the partitioned state ordering
    P = ml.zeros((N,N))
    for i in range(Nz):
        P[i,ixz[i]]=1
    for i in range(Np):
        P[Nz+i,ixp[i]]=1
    for i in range(Nn):
        P[Nz+Np+i,ixn[i]]=1
    iP = P.I
    Qv = P*Q*iP
    Rv = P*R*iP

    # new fluid process censored to states + and -
    iQv00 = la.pinv(-Qv[:Nz,:Nz])
    Qbar = Qv[Nz:, Nz:] + Qv[Nz:,:Nz]*iQv00*Qv[:Nz,Nz:]
    absRi = Diag(np.abs(1./np.diag(Rv[Nz:,Nz:])))
    Qz = absRi * Qbar

    Psi, K, U = FluidFundamentalMatrices (Qz[:Np,:Np], Qz[:Np,Np:], Qz[Np:,:Np], Qz[Np:,Np:], "PKU", prec)

    # closing matrix
    Pm = np.hstack((ml.eye(Np), Psi)) * absRi
    iCn = absRi[Np:,Np:]
    iCp = absRi[:Np,:Np]
    clo = np.hstack(((iCp*Qv[Nz:Nz+Np,:Nz]+Psi*iCn*Qv[Nz+Np:,:Nz])*iQv00, Pm))
    
    if len(Q0)==0: # regular boundary behavior
        clo = clo * P # go back the the original state ordering

        # calculate boundary vector   
        Ua = iCn*Qv[Nz+Np:,:Nz]*iQv00*ml.ones((Nz,1)) + iCn*ml.ones((Nn,1)) + Qz[Np:,:Np]*la.inv(-K)*clo*ml.ones((Nz+Np+Nn,1))
        pm = Linsolve (ml.hstack((U,Ua)).T, ml.hstack((ml.zeros((1,Nn)),ml.ones((1,1)))).T).T

        # create the result
        mass0 = ml.hstack((pm*iCn*Qv[Nz+Np:,:Nz]*iQv00, ml.zeros((1,Np)), pm*iCn))*P
        ini = pm*Qz[Np:,:Np]        
    else:
        # solve a linear system for ini(+), pm(-) and pm(0)        
        Q0v = P*Q0*iP
        M = ml.vstack((-clo*Rv, Q0v[Nz+Np:,:], Q0v[:Nz,:]))
        Ma = ml.vstack((np.sum(la.inv(-K)*clo,1), ml.ones((Nz+Nn,1))))
        sol = Linsolve (ml.hstack((M,Ma)).T, ml.hstack((ml.zeros((1,N)),ml.ones((1,1)))).T).T;
        ini = sol[:,:Np]
        clo = clo * P
        mass0 = ml.hstack((sol[:,Np+Nn:], ml.zeros((1,Np)), sol[:,Np:Np+Nn]))*P

    return mass0, ini, K, clo

def FluidStationaryDistr (mass0, ini, K, clo, x):
    """
    Returns the stationary distribution of a Markovian 
    fluid model at the given points.
    
    Parameters
    ----------
    mass0 : matrix, shape (1,Np+Nm)
        The stationary probability vector of zero level
    ini : matrix, shape (1,Np)
        The initial vector of the stationary density
    K : matrix, shape (Np,Np)
        The matrix parameter of the stationary density
    clo : matrix, shape (Np,Np+Nm)
        The closing matrix of the stationary density
    x : vector, length (K)
        The distribution function is computed at these 
        points.
    
    Returns
    -------
    pi : matrix, shape (K,Nm+Np)
        The ith row of pi is the probability that the fluid
        level is less than or equal to x(i), while being in
        different states of the background process.
    """

    m = clo.shape[1]
    y = ml.empty((len(x),m))
    closing = -K.I*clo
    for i in range(len(x)):
        y[i,:] = mass0 + ini*(ml.eye(K.shape[0])-la.expm(K*x[i]))*closing

    return y
