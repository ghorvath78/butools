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

def MarginalDistributionFromDRAP (H0, H1):
    """
    Returns the matrix geometrically distributed marginal 
    distribution of a discrete rational arrival process.
    
    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the discrete rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the discrete rational arrival process
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14
    
    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial vector of the matrix geometrically
        distributed marginal distribution
    A : matrix, shape (M,M)
        The matrix parameter of the matrix geometrically
        distributed marginal distribution    
    """

    if butools.checkInput and not CheckDRAPRepresentation (H0, H1):
        raise Exception("MarginalDistributionFromDRAP: Input is not a valid DRAP representation!")    

    return (DRPSolve(la.inv(ml.eye(H0.shape[0])-H0)*H1), H0)

def MarginalDistributionFromDMAP (D0, D1):
    """
    Returns the discrete phase type distributed marginal 
    distribution of a discrete Markovian arrival process.
    
    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the discrete Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the discrete Markovian arrival process
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14
    
    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial probability vector of the discrete phase 
        type distributed marginal distribution
    A : matrix, shape (M,M)
        The transient generator of the discrete phase type 
        distributed marginal distribution    
    """

    if butools.checkInput and not CheckDMAPRepresentation (D0, D1):
        raise Exception("MarginalDistributionFromDMAP: Input is not a valid DMAP representation!")    

    return MarginalDistributionFromDRAP (D0, D1)

def MarginalDistributionFromDMRAP (H):
    """
    Returns the matrix geometrically distributed marginal 
    distribution of a discrete marked rational arrival 
    process.
    
    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the DMRAP
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14
    
    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial vector of the matrix geometrically
        distributed marginal distribution
    A : matrix, shape (M,M)
        The matrix parameter of the matrix geometrically
        distributed marginal distribution    
    """

    if butools.checkInput and not CheckDMRAPRepresentation (H):
        raise Exception("MarginalDistributionFromDMRAP: Input is not a valid DMRAP representation!")    

    return (DRPSolve(la.inv(ml.eye(H[0].shape[0])-H[0])*SumMatrixList(H[1:])), H[0])

def MarginalDistributionFromDMMAP (D):
    """
    Returns the discrete phase type distributed marginal 
    distribution of a discrete marked Markovian arrival 
    process.
    
    Parameters
    ----------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the DMMAP
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14
    
    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial probability vector of the discrete phase 
        type distributed marginal distribution
    A : matrix, shape (M,M)
        The transient generator of the discrete phase type 
        distributed marginal distribution    
    """

    if butools.checkInput and not CheckDMMAPRepresentation (D):
        raise Exception("MarginalDistributionFromDMMAP: Input is not a valid DMMAP representation!")    

    return MarginalDistributionFromDMRAP (D)

from butools.dph import MomentsFromDPH, MomentsFromMG

def MarginalMomentsFromDRAP (H0, H1, K=0):
    """
    Returns the moments of the marginal distribution of a 
    discrete rational arrival process.
    
    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the discrete rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the discrete rational arrival process
    K : int, optional
        Number of moments to compute. If K=0, 2*M-1 moments
        are computed. The default value is K=0.
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14
    
    Returns
    -------
    moms : row vector of doubles, length K
        The vector of moments.
    """

    if butools.checkInput and not CheckDRAPRepresentation (H0, H1):
        raise Exception("MarginalMomentsFromDRAP: Input is not a valid DRAP representation!")    

    alpha,A = MarginalDistributionFromDRAP(H0,H1)
    return MomentsFromMG(alpha,A,K)

def MarginalMomentsFromDMAP (D0, D1, K=0):
    """
    Returns the moments of the marginal distribution of a 
    discrete Markovian arrival process.
    
    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the discrete Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the discrete Markovian arrival process
    K : int, optional
        Number of moments to compute. If K=0, 2*M-1 moments
        are computed. The default value is K=0.
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14
    
    Returns
    -------
    moms : row vector of doubles, length K
        The vector of moments.
    """

    if butools.checkInput and not CheckDMAPRepresentation (D0, D1):
        raise Exception("MarginalMomentsFromDMAP: Input is not a valid DMAP representation!")    

    alpha,A = MarginalDistributionFromDMAP(D0,D1)
    return MomentsFromDPH(alpha,A,K)
    
def MarginalMomentsFromDMRAP (H, K=0):
    """
    Returns the moments of the marginal distribution of a 
    discrete marked rational arrival process.
    
    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the DMRAP
    K : int, optional
        Number of moments to compute. If K=0, 2*M-1 moments
        are computed. The default value is K=0.
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14
    
    Returns
    -------
    moms : row vector of doubles, length K
        The vector of moments.
    """

    if butools.checkInput and not CheckDMRAPRepresentation (H):
        raise Exception("MarginalMomentsFromDMRAP: Input is not a valid DMRAP representation!")    

    alpha,A = MarginalDistributionFromDMRAP(H)
    return MomentsFromMG(alpha,A,K)

def MarginalMomentsFromDMMAP (D, K=0):
    """
    Returns the moments of the marginal distribution of a 
    discrete marked Markovian arrival process.
    
    Parameters
    ----------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the DMMAP
    K : int, optional
        Number of moments to compute. If K=0, 2*M-1 moments
        are computed. The default value is K=0.
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14
    
    Returns
    -------
    moms : row vector of doubles, length K
        The vector of moments.
    """

    if butools.checkInput and not CheckDMMAPRepresentation (D):
        raise Exception("MarginalMomentsFromDMMAP: Input is not a valid DMMAP representation!")    

    alpha,A = MarginalDistributionFromDMMAP(D)
    return MomentsFromDPH(alpha,A,K)

def LagCorrelationsFromDRAP (H0, H1, L=1):
    """
    Returns the lag autocorrelations of a discrete rational 
    arrival process.
    
    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the discrete rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the discrete rational arrival process
    L : double, optional
        The number of lags to compute. The default value is 1
    prec : double, optional
        Numerical precision to check if the input is valid. 
        The default value is 1e-14
    
    Returns
    -------
    acf : column vector of doubles, length (L)
        The lag autocorrelation function up to lag L
        
    """

    if butools.checkInput and not CheckDRAPRepresentation (H0, H1):
        raise Exception("LagCorrelationsFromDRAP: Input is not a valid DRAP representation!")    

    H0i = la.inv(ml.eye(H0.shape[0])-H0)
    P = H0i*H1
    pi = DRPSolve(P)
    m1, m2 = MomentsFromMG(pi, H0, 2)
    pi = pi * H0i * P

    corr = []
    for i in range(L):
        corr.append((np.sum(pi*H0i) - m1*m1) / (m2 - m1*m1))
        pi = pi * P
#    if L>1:
    return np.array(corr)
#    else:
#        return corr[0]

def LagCorrelationsFromDMAP (D0, D1, L=1):
    """
    Returns the lag autocorrelations of a discrete Markovian 
    arrival process.
    
    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the discrete Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the discrete Markovian arrival process
    L : double, optional
        The number of lags to compute. The default value is 1
    prec : double, optional
        Numerical precision to check if the input is valid. 
        The default value is 1e-14
    
    Returns
    -------
    acf : column vector of doubles, length (L)
        The lag autocorrelation function up to lag L
        
    """

    if butools.checkInput and not CheckDMAPRepresentation (D0, D1):
        raise Exception("LagCorrelationsFromDMAP: Input is not a valid DMAP representation!")    

    return LagCorrelationsFromDRAP (D0, D1, L)

def LagkJointMomentsFromDMRAP (H, K=0, L=1):
    """
    Returns the lag-L joint moments of a discrete marked 
    rational arrival process.
    
    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the DMRAP to check
    K : int, optional
        The dimension of the matrix of joint moments to 
        compute. If K=0, the MxM joint moments will be 
        computed. The default value is 0
    L : int, optional
        The lag at which the joint moments are computed.
        The default value is 1
    prec : double, optional
        Numerical precision to check if the input is valid. 
        The default value is 1e-14
    
    Returns
    -------
    Nm : list/cell of matrices of shape(K+1,K+1), length(L)
        The matrices containing the lag-L joint moments,
        starting from moment 0.
    """

    if butools.checkInput and not CheckDMRAPRepresentation (H):
        raise Exception("LagkJointMomentsFromDMRAP: Input is not a valid DMRAP representation!")    

    if K==0:
        K = H[0].shape[0]-1
    M = len(H)-1
    H0 = H[0]
    sumH = SumMatrixList(H[1:])

    H0i = la.inv(ml.eye(H0.shape[0])-H0)
    P = H0i*sumH
    pi = DRPSolve(P)
    
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

def LagkJointMomentsFromDMMAP (D, K=0, L=1):
    """
    Returns the lag-L joint moments of a discrete marked 
    Markovian arrival process.
    
    Parameters
    ----------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the DMMAP to check
    K : int, optional
        The dimension of the matrix of joint moments to 
        compute. If K=0, the MxM joint moments will be 
        computed. The default value is 0
    L : int, optional
        The lag at which the joint moments are computed.
        The default value is 1
    prec : double, optional
        Numerical precision to check if the input is valid. 
        The default value is 1e-14
    
    Returns
    -------
    Nm : list/cell of matrices of shape(K+1,K+1), length(L)
        The matrices containing the lag-L joint moments,
        starting from moment 0.
    """

    if butools.checkInput and not CheckDMMAPRepresentation (D):
        raise Exception("LagkJointMomentsFromDMMAP: Input is not a valid DMMAP representation!")    

    return LagkJointMomentsFromDMRAP(D, K, L)

    
def LagkJointMomentsFromDRAP (H0, H1, K=0, L=1):
    """
    Returns the lag-L joint moments of a discrete rational 
    arrival process.
    
    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the discrete rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the discrete rational arrival process
    K : int, optional
        The dimension of the matrix of joint moments to 
        compute. If K=0, the MxM joint moments will be 
        computed. The default value is 0
    L : int, optional
        The lag at which the joint moments are computed.
        The default value is 1
    prec : double, optional
        Numerical precision to check if the input is valid. 
        The default value is 1e-14
    
    Returns
    -------
    Nm : matrix, shape(K+1,K+1)
        The matrix containing the lag-L joint moments, 
        starting from moment 0.
    """

    if butools.checkInput and not CheckDRAPRepresentation (H0, H1):
        raise Exception("LagkJointMomentsFromDRAP: Input is not a valid DRAP representation!")    

    return LagkJointMomentsFromDMRAP((H0,H1), K, L)[0]

def LagkJointMomentsFromDMAP (D0, D1, K=0, L=1):
    """
    Returns the lag-L joint moments of a discrete Markovian 
    arrival process.
    
    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the discrete Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the discrete Markovian arrival process
    K : int, optional
        The dimension of the matrix of joint moments to 
        compute. If K=0, the MxM joint moments will be 
        computed. The default value is 0
    L : int, optional
        The lag at which the joint moments are computed.
        The default value is 1
    prec : double, optional
        Numerical precision to check if the input is valid. 
        The default value is 1e-14
    
    Returns
    -------
    Nm : matrix, shape(K+1,K+1)
        The matrix containing the lag-L joint moments, 
        starting from moment 0.
    """

    if butools.checkInput and not CheckDMAPRepresentation (D0, D1):
        raise Exception("LagkJointMomentsFromDMAP: Input is not a valid DMAP representation!")    

    return LagkJointMomentsFromDRAP(D0, D1, K, L)
    
def RandomDMMAP (order, types, mean=10.0, zeroEntries=0, maxTrials=1000, prec=1e-7):
    """
    Returns a random discrete Markovian arrival process.
    
    Parameters
    ----------
    order : int
        The size of the DMAP
    mean : double, optional
        The mean inter-arrival times of the DMMAP
    types : int
        The number of different arrival types
    zeroEntries : int, optional
        The number of zero entries in the D0 and D1 matrices
    maxTrials : int, optional
        The maximum number of trials to find a proper DMMAP 
        (that has an irreducible phase process and none of 
        its parameters is all-zero)
    prec : double, optional
        Numerical precision for checking the irreducibility.
        The default value is 1e-14.
    
    Returns
    -------
    D : list/cell of matrices of shape(M,M), length(types+1)
        The D0...Dtypes matrices of the DMMAP 
    
    Notes
    -----
    If it fails, try to increase the 'maxTrials' parameter,
    or/and the 'mean' parameter.
    """

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
                if np.min(np.abs(alpha)) > prec:
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
                            m = MarginalMomentsFromDMMAP(Dv, 1)[0]
                            d = 1 - (1-d)*m/mean
                            for i in range(types+1):         
                                D[i] = Diag(1-d)*D[i]
                            D[0] = D[0] + Diag(d)                        
                            if CheckDMMAPRepresentation(D):
                                return D
                        except:
                            pass
            trials += 1
    raise Exception("No feasible random DMAP/DMMAP found with such many zero entries! Try to increase the maxTrials parameter!")

def RandomDMAP (order, mean=10.0, zeroEntries=0, maxTrials=1000, prec=1e-7):
    """
    Returns a random disctere Markovian arrival process.
    
    Parameters
    ----------
    order : int
        The size of the DMAP
    mean : double, optional
        The mean inter-arrival times of the DMAP
    zeroEntries : int, optional
        The number of zero entries in the D0 and D1 matrices
    maxTrials : int, optional
        The maximum number of trials to find a proper DMAP 
        (that has an irreducible phase process and none of 
        its parameters is all-zero)
    prec : double, optional
        Numerical precision for checking the irreducibility.
        The default value is 1e-14.
    
    Returns
    -------
    D0 : vector, shape (1,M)
        The D0 matrix of the DMAP
    D1 : matrix, shape (M,M)
        The D1 matrix of the DMAP
    
    Notes
    -----
    If it fails, try to increase the 'maxTrials' parameter,
    or/and the 'mean' parameter.
    """
    
    return RandomDMMAP (order, 1, mean, zeroEntries, maxTrials, prec)
