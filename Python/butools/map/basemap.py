# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 10:42:38 2013

@author: gabor
"""

from math import sqrt
import numpy as np
import numpy.matlib as ml
import numpy.linalg as la
import butools
from butools.utils import Diag
from butools.mc import CTMCSolve, DRPSolve
from butools.map import CheckRAPRepresentation, CheckMAPRepresentation, CheckMRAPRepresentation, CheckMMAPRepresentation

def MarginalDistributionFromRAP (H0, H1, prec=1e-14):

    if butools.checkInput and not CheckRAPRepresentation (H0, H1, prec):
        raise Exception("MarginalDistributionFromRAP: Input is not a valid RAP representation!")    

    return (DRPSolve(la.inv(-H0)*H1, prec), H0)

def MarginalDistributionFromMAP (D0, D1, prec=1e-14):

    if butools.checkInput and not CheckMAPRepresentation (D0, D1, prec):
        raise Exception("MarginalDistributionFromMAP: Input is not a valid MAP representation!")    

    return MarginalDistributionFromRAP (D0, D1, prec)

def MarginalDistributionFromMRAP (H, prec=1e-14):

    if butools.checkInput and not CheckMRAPRepresentation (H, prec):
        raise Exception("MarginalDistributionFromMRAP: Input is not a valid MRAP representation!")    

    Hk = ml.matrix(H[1])
    for i in range (2, len(H)):
        Hk += H[i]
    return (DRPSolve(la.inv(-H[0])*Hk, prec), H[0])

def MarginalDistributionFromMMAP (D, prec=1e-14):

    if butools.checkInput and not CheckMMAPRepresentation (D, prec):
        raise Exception("MarginalDistributionFromMMAP: Input is not a valid MMAP representation!")    

    return MarginalDistributionFromMRAP (D, prec)

from butools.ph import MomentsFromPH, MomentsFromME

def MarginalMomentsFromRAP (H0, H1, K=0, prec=1e-14):

    if butools.checkInput and not CheckRAPRepresentation (H0, H1, prec):
        raise Exception("MarginalMomentsFromRAP: Input is not a valid RAP representation!")    

    alpha,A = MarginalDistributionFromRAP(H0,H1,prec)
    return MomentsFromME(alpha,A,K,prec)

def MarginalMomentsFromMAP (D0, D1, K=0, prec=1e-14):

    if butools.checkInput and not CheckMAPRepresentation (D0, D1, prec):
        raise Exception("MarginalMomentsFromMAP: Input is not a valid MAP representation!")    

    alpha,A = MarginalDistributionFromMAP(D0,D1,prec)
    return MomentsFromPH(alpha,A,K,prec)
    
def MarginalMomentsFromMRAP (H, K=0, prec=1e-14):

    if butools.checkInput and not CheckMRAPRepresentation (H, prec):
        raise Exception("MarginalMomentsFromMRAP: Input is not a valid MRAP representation!")    

    alpha,A = MarginalDistributionFromMRAP(H,prec)
    return MomentsFromME(alpha,A,K,prec)

def MarginalMomentsFromMMAP (D, K=0, prec=1e-14):

    if butools.checkInput and not CheckMMAPRepresentation (D, prec):
        raise Exception("MarginalMomentsFromMMAP: Input is not a valid MMAP representation!")    

    alpha,A = MarginalDistributionFromMMAP(D,prec)
    return MomentsFromPH(alpha,A,K,prec)

def LagCorrelationsFromRAP (H0, H1, L=1, prec=1e-14):

    if butools.checkInput and not CheckRAPRepresentation (H0, H1, prec):
        raise Exception("LagCorrelationsFromRAP: Input is not a valid RAP representation!")    

    H0i = la.inv(-H0)
    P = H0i*H1
    pi = DRPSolve(P, prec)
    m1, m2 = MomentsFromME(pi, H0, 2, prec)
    pi = pi * H0i * P

    corr = []
    for i in range(L):
        corr.append((np.sum(pi*H0i) - m1*m1) / (m2 - m1*m1))
        pi = pi * P
    if L>1:
        return corr
    else:
        return corr[0]

def LagCorrelationsFromMAP (D0, D1, L=1, prec=1e-14):

    if butools.checkInput and not CheckMAPRepresentation (D0, D1, prec):
        raise Exception("LagCorrelationsFromMAP: Input is not a valid MAP representation!")    

    return LagCorrelationsFromRAP (D0, D1, L, prec)

def LagkJointMomentsFromMRAP (H, K=0, L=1, prec=1e-14):

    if butools.checkInput and not CheckMRAPRepresentation (H, prec):
        raise Exception("LagkJointMomentsFromMRAP: Input is not a valid MRAP representation!")    

    if K==0:
        K = H[0].shape[0]-1
    M = len(H)-1
    H0 = H[0]
    sumH = ml.zeros(H[0].shape)
    for i in range(M):
        sumH += H[i+1]

    H0i = la.inv(-H0)
    P = H0i*sumH
    pi = DRPSolve(P, prec)
    
    Pw = ml.eye(H0.shape[0])
    H0p = [ml.matrix(Pw)]
    for i in range(1,K+1):
        Pw *= i*H0i
        H0p.append(ml.matrix(Pw))

    Pl = la.matrix_power (P, L-1)

    Nm = []
    for m in range(M):
        Nmm = ml.zeros ((K+1,K+1))
        for i in range(K+1):
            for j in range(K+1):
                Nmm[i,j] = np.sum (pi * H0p[i] * H0i * H[m+1] * Pl * H0p[j])
        Nm.append(ml.matrix(Nmm))
    return Nm

def LagkJointMomentsFromMMAP (D, K=0, L=1, prec=1e-14):

    if butools.checkInput and not CheckMMAPRepresentation (D, prec):
        raise Exception("LagkJointMomentsFromMMAP: Input is not a valid MMAP representation!")    

    return LagkJointMomentsFromMRAP(D, K, L, prec)

    
def LagkJointMomentsFromRAP (H0, H1, K=0, L=1, prec=1e-14):

    if butools.checkInput and not CheckRAPRepresentation (H0, H1, prec):
        raise Exception("LagkJointMomentsFromRAP: Input is not a valid RAP representation!")    

    return LagkJointMomentsFromMRAP((H0,H1), K, L, prec)[0]

def LagkJointMomentsFromMAP (D0, D1, K=0, L=1, prec=1e-14):

    if butools.checkInput and not CheckMAPRepresentation (D0, D1, prec):
        raise Exception("LagkJointMomentsFromMAP: Input is not a valid MAP representation!")    

    return LagkJointMomentsFromRAP(D0, D1, K, L, prec)
    
def RandomMMAP (order, types, mean=1.0, zeroEntries=0, maxTrials=1000, prec=1e-14):
    # distribute the zero entries among the rows
    def allZeroDistr (states, zeros):
        if states==1:
            return [[zeros]]
        else:
            o = [];
            for i in range(zeros+1):
                x = allZeroDistr (states-1, zeros-i)
                for j in range(len(x)):
                    xt = x[j]
                    xt.append(i)
                    xt.sort()
                    # check if we have it already
                    if o.count(xt)==0:
                        o.append(xt)
            return o

    if zeroEntries > (types+1)*order*order - 2*order:
        raise Exception("RandomMAP/MMAP: Too many zero entries requested! Try to decrease the zeroEntries parameter!")

    zeroDistr = allZeroDistr(order, zeroEntries)   

    trials = 1
    while trials<maxTrials:
        # select a configuration from zeroDistr: it is a list describing the zero entries in each row
        zdix = np.random.permutation(len(zeroDistr))
        for k in range(len(zeroDistr)):
            zDistr = zeroDistr[zdix[k]];
            bad = False            
            for d in zDistr:
                if d>=(types+1)*order-1:
                    bad = True
                    break
            if bad:
                continue
            B = np.zeros((order,(types+1)*order))
            for i in range(order):
                rp = np.random.permutation((types+1)*order-1)
                a = np.zeros((types+1)*order-1)
                for j in range((types+1)*order - 1 - zDistr[i]):
                    a[rp[j]] = np.random.rand()
                B[i,0:i] = a[0:i]
                B[i,i+1:] = a[i:]
            # construct MMAP matrices
            D = []
            for i in range(types+1):         
                Di = ml.matrix(B[:,i*order:(i+1)*order])
                D.append(Di)           
            D[0] -= Diag(np.sum(B,1))
            # check if it is a proper MAP (irreducible phase process & no full zero matrix)
            sumD = ml.zeros((order,order))
            for i in range(len(D)):         
                sumD += D[i]
            if la.matrix_rank(D[0])==order and la.matrix_rank(sumD)==order-1:
                alpha = CTMCSolve(sumD)
                if np.min(np.abs(alpha)) > sqrt(prec):
                    fullZero = False
                    for Di in D:
                        if np.all(Di==0.0):
                            fullZero = True
                            break
                    if not fullZero:
                        # scale to the mean value
                        m = MarginalMomentsFromMMAP (D, 1, prec)[0]
                        for i in range(types+1):
                            D[i] *= m / mean
                        return D
            trials += 1
    raise Exception("No feasible random MAP/MMAP found with such many zero entries! Try to increase the maxTrials parameter!")

def RandomMAP (order, mean=1.0, zeroEntries=0, maxTrials=1000, prec=1e-14):
    
    return RandomMMAP (order, 1, mean, zeroEntries, maxTrials, prec)
