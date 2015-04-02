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
    
    Psi, K, U = FluidFundamentalMatrices (Fpp, Fpm, Fmp, Fmm, "PKU", prec)
    
    mass0 = CTMCSolve(U, prec)
    nr = np.sum(mass0) + 2.0*np.sum(mass0*Fmp*-K.I)
    
    return mass0/nr, mass0*Fmp, K, ml.hstack((ml.eye(Fpp.shape[0]), Psi))

def GeneralFluidSolve (Q, R, Q0=None, prec=1e-14):
    
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
    iQv00 = np.pinv(-Qv[:Nz,:Nz])
    Qbar = Qv[Nz:, Nz:] + Qv[Nz:,:Nz]*iQv00*Qv[:Nz,Nz:]
    absRi = Diag(np.abs(1./np.diag(Rv[Nz:,Nz:])))
    Qz = absRi * Qbar

    Psi, K, U = FluidFundamentalMatrices (Qz[:Np,:Np], Qz[:Np,Np:], Qz[Np:,:Np], Qz[Np:,Np:], "PKU")

    # closing matrix
    Pm = np.hstack((ml.eye(Np), Psi)) * absRi
    iCn = absRi[Np:,Np:]
    iCp = absRi[:Np,:Np]
    clo = np.hstack(((iCp*Qv[Nz:,Nz+Np,:Nz]+Psi*iCn*Qv[Nz+Np:,:Nz])*iQv00, Pm))
    
    if Q0==None: # regular boundary behavior
        clo = clo * P # go back the the original state ordering

        # calculate boundary vector   
        Ua = iCn*Qv[Nz+Np:,:Nz]*iQv00*ml.ones((Nz,1)) + iCn*ml.ones((Nn,1)) + Qz[Np:,:Np]*la.inv(-K)*clo*ml.ones((Nz+Np+Nn,1))
        pm = Linsolve (ml.hstack((U,Ua)).T, ml.hstack((ml.zeros((1,Nn)),1)).T).T

        # create the result
        mass0 = ml.hstack((pm*iCn*Qv[Nz+Np:,:Nz]*iQv00, ml.zeros((1,Np)), pm*iCn))*P
        ini = pm*Qz[Np:,:Np]        
    else:
        # solve a linear system for ini(+), pm(-) and pm(0)        
        Q0v = P*Q0*iP
        M = ml.vstack((-clo*Rv, Q0v[Nz+Np:,:], Q0v[:Nz,:]))
        Ma = ml.vstack((np.sum(la.inv(-K)*clo,1), ml.ones((Nz+Nn,1))))
        sol = Linsolve (ml.hstack((M,Ma)).T, ml.hstack((ml.zeros((1,N)),1)).T).T;
        ini = sol[:Np]
        clo = clo * P
        mass0 = ml.hstack((sol[Np+Nn:], ml.zeros((1,Np)), sol[Np:Np+Nn]))*P

    return mass0, ini, K, clo

def FluidStationaryDistr (mass0, ini, K, clo, x):

    m = clo.shape[1]
    y = ml.empty((len(x),m))
    closing = -K.I*clo
    for i in range(len(x)):
        y[i,:] = mass0 + ini*(ml.eye(K.shape[0])-la.expm(K*x[i]))*closing

    return y
