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
from butools.utils import Diag, SumMatrixList
from butools.moments import MomsFromFactorialMoms, JMomsFromJFactorialMoms
from butools.mc import DTMCSolve, DRPSolve
from butools.dmap import CheckDRAPRepresentation, CheckDMAPRepresentation, CheckDMRAPRepresentation, CheckDMMAPRepresentation

def MarginalDistributionFromDRAP (H0, H1, prec=1e-14):

    if butools.checkInput and not CheckDRAPRepresentation (H0, H1, prec):
        raise Exception("MarginalDistributionFromDRAP: Input is not a valid DRAP representation!")    

    return (DRPSolve(la.inv(ml.eye(H0.shape[0])-H0)*H1, prec), H0)

def MarginalDistributionFromDMAP (D0, D1, prec=1e-14):

    if butools.checkInput and not CheckDMAPRepresentation (D0, D1, prec):
        raise Exception("MarginalDistributionFromDMAP: Input is not a valid DMAP representation!")    

    return MarginalDistributionFromDRAP (D0, D1, prec)

def MarginalDistributionFromDMRAP (H, prec=1e-14):

    if butools.checkInput and not CheckDMRAPRepresentation (H, prec):
        raise Exception("MarginalDistributionFromDMRAP: Input is not a valid DMRAP representation!")    

    return (DRPSolve(la.inv(ml.eye(H[0].shape[0])-H[0])*SumMatrixList(H[1:]), prec), H[0])

def MarginalDistributionFromDMMAP (D, prec=1e-14):

    if butools.checkInput and not CheckDMMAPRepresentation (D, prec):
        raise Exception("MarginalDistributionFromDMMAP: Input is not a valid DMMAP representation!")    

    return MarginalDistributionFromDMRAP (D, prec)

from butools.dph import MomentsFromDPH, MomentsFromMG

def MarginalMomentsFromDRAP (H0, H1, K=0, prec=1e-14):

    if butools.checkInput and not CheckDRAPRepresentation (H0, H1, prec):
        raise Exception("MarginalMomentsFromDRAP: Input is not a valid DRAP representation!")    

    alpha,A = MarginalDistributionFromDRAP(H0,H1,prec)
    return MomentsFromMG(alpha,A,K,prec)

def MarginalMomentsFromDMAP (D0, D1, K=0, prec=1e-14):

    if butools.checkInput and not CheckDMAPRepresentation (D0, D1, prec):
        raise Exception("MarginalMomentsFromDMAP: Input is not a valid DMAP representation!")    

    alpha,A = MarginalDistributionFromDMAP(D0,D1,prec)
    return MomentsFromDPH(alpha,A,K,prec)
    
def MarginalMomentsFromDMRAP (H, K=0, prec=1e-14):

    if butools.checkInput and not CheckDMRAPRepresentation (H, prec):
        raise Exception("MarginalMomentsFromDMRAP: Input is not a valid DMRAP representation!")    

    alpha,A = MarginalDistributionFromDMRAP(H,prec)
    return MomentsFromMG(alpha,A,K,prec)

def MarginalMomentsFromDMMAP (D, K=0, prec=1e-14):

    if butools.checkInput and not CheckDMMAPRepresentation (D, prec):
        raise Exception("MarginalMomentsFromDMMAP: Input is not a valid DMMAP representation!")    

    alpha,A = MarginalDistributionFromDMMAP(D,prec)
    return MomentsFromDPH(alpha,A,K,prec)

def LagCorrelationsFromDRAP (H0, H1, L=1, prec=1e-14):

    if butools.checkInput and not CheckDRAPRepresentation (H0, H1, prec):
        raise Exception("LagCorrelationsFromDRAP: Input is not a valid DRAP representation!")    

    H0i = la.inv(ml.eye(H0.shape[0])-H0)
    P = H0i*H1
    pi = DRPSolve(P, prec)
    m1, m2 = MomentsFromMG(pi, H0, 2, prec)
    pi = pi * H0i * P

    corr = []
    for i in range(L):
        corr.append((np.sum(pi*H0i) - m1*m1) / (m2 - m1*m1))
        pi = pi * P
    if L>1:
        return corr
    else:
        return corr[0]

def LagCorrelationsFromDMAP (D0, D1, L=1, prec=1e-14):

    if butools.checkInput and not CheckDMAPRepresentation (D0, D1, prec):
        raise Exception("LagCorrelationsFromDMAP: Input is not a valid DMAP representation!")    

    return LagCorrelationsFromDRAP (D0, D1, L, prec)

def LagkJointMomentsFromDMRAP (H, K=0, L=1, prec=1e-14):

    if butools.checkInput and not CheckDMRAPRepresentation (H, prec):
        raise Exception("LagkJointMomentsFromDMRAP: Input is not a valid DMRAP representation!")    

    if K==0:
        K = H[0].shape[0]-1
    M = len(H)-1
    H0 = H[0]
    sumH = SumMatrixList(H[1:])

    H0i = la.inv(ml.eye(H0.shape[0])-H0)
    P = H0i*sumH
    pi = DRPSolve(P, prec)
    
    Pw = ml.eye(H0.shape[0])
    H0p = [ml.matrix(Pw)]
    Pw = Pw*H0i
    H0p.append(Pw)
    for i in range(2,K+1):
        Pw = Pw*i*H0i*H0
        H0p.append(ml.matrix(Pw))

    Pl = la.matrix_power (P, L-1)

    Nm = []
    for m in range(M):
        Nmm = ml.zeros ((K+1,K+1))
        for i in range(K+1):
            for j in range(K+1):
                Nmm[i,j] = np.sum (pi * H0p[i] * H0i * H[m+1] * Pl * H0p[j])
        row1 = np.array([MomsFromFactorialMoms(Nmm[0,1:].A.flatten())])
        col1 = np.array([MomsFromFactorialMoms(Nmm[1:,0].A.flatten())]).T
        mid = JMomsFromJFactorialMoms(Nmm[1:,1:])
        Nm.append(np.bmat([[[[Nmm[0,0]]], row1 ], [col1, mid]]))
    return Nm

def LagkJointMomentsFromDMMAP (D, K=0, L=1, prec=1e-14):

    if butools.checkInput and not CheckDMMAPRepresentation (D, prec):
        raise Exception("LagkJointMomentsFromDMMAP: Input is not a valid DMMAP representation!")    

    return LagkJointMomentsFromDMRAP(D, K, L, prec)

    
def LagkJointMomentsFromDRAP (H0, H1, K=0, L=1, prec=1e-14):

    if butools.checkInput and not CheckDRAPRepresentation (H0, H1, prec):
        raise Exception("LagkJointMomentsFromDRAP: Input is not a valid DRAP representation!")    

    return LagkJointMomentsFromDMRAP((H0,H1), K, L, prec)[0]

def LagkJointMomentsFromDMAP (D0, D1, K=0, L=1, prec=1e-14):

    if butools.checkInput and not CheckDMAPRepresentation (D0, D1, prec):
        raise Exception("LagkJointMomentsFromDMAP: Input is not a valid DMAP representation!")    

    return LagkJointMomentsFromDRAP(D0, D1, K, L, prec)
    
def RandomDMMAP (order, types, mean=1.0, zeroEntries=0, maxTrials=1000, prec=1e-14):

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
        raise Exception("RandomDMAP/DMMAP: Too many zero entries requested! Try to decrease the zeroEntries parameter!")

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
            # construct DMMAP matrices
            D = []
            sc = np.zeros(order)
            for i in range(types+1):         
                Di = ml.matrix(B[:,i*order:(i+1)*order])
                D.append(Di)
                sc += np.sum(Di,1).A.flatten()
            if np.any(sc==0):
                continue
            for i in range(types+1):
                D[i] = Diag(1.0/sc) * D[i]
            # check if it is a proper DMAP (irreducible phase process & no full zero matrix)
            sumD = SumMatrixList(D)
            if la.matrix_rank(D[0])==order and la.matrix_rank(ml.eye(D[0].shape[0])-sumD)==order-1:
                alpha = DTMCSolve(sumD)
                if np.min(np.abs(alpha)) > sqrt(prec):
                    fullZero = False
                    for Di in D:
                        if np.all(Di==0.0):
                            fullZero = True
                            break
                    if not fullZero:
                        # diagonals of matrix A:
                        d = np.random.rand(order)
                        # scale to the mean value
                        Dv = []
                        for i in range(types+1):         
                            Dv.append(Diag(1-d)*D[i])
                        Dv[0] = Dv[0] + Diag(d)                        
                        try:                        
                            m = MarginalMomentsFromDMMAP(Dv, 1, prec)[0]
                            d = 1 - (1-d)*m/mean
                            for i in range(types+1):         
                                D[i] = Diag(1-d)*D[i]
                            D[0] = D[0] + Diag(d)                        
                            if CheckDMMAPRepresentation(D,prec):
                                return D
                        except:
                            pass
            trials += 1
    raise Exception("No feasible random DMAP/DMMAP found with such many zero entries! Try to increase the maxTrials parameter!")

def RandomDMAP (order, mean=1.0, zeroEntries=0, maxTrials=1000, prec=1e-14):
    
    return RandomDMMAP (order, 1, mean, zeroEntries, maxTrials, prec)
