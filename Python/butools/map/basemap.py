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

def MarginalDistributionFromRAP (H0, H1):
    """
    Returns the phase type distributed marginal distribution
    of a rational arrival process.
    
    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational arrival process
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14
    
    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial vector of the matrix exponentially
        distributed marginal
    A : matrix, shape (M,M)
        The matrix parameter of the matrix exponentially
        distributed marginal    
    """

    if butools.checkInput and not CheckRAPRepresentation (H0, H1):
        raise Exception("MarginalDistributionFromRAP: Input is not a valid RAP representation!")    

    return (DRPSolve(la.inv(-H0)*H1), H0)

def MarginalDistributionFromMAP (D0, D1):
    """
    Returns the phase type distributed marginal distribution
    of a Markovian arrival process.
    
    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14
    
    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial probability vector of the phase type 
        distributed marginal distribution
    A : matrix, shape (M,M)
        The transient generator of the phase type distributed
        marginal distribution    
    """

    if butools.checkInput and not CheckMAPRepresentation (D0, D1):
        raise Exception("MarginalDistributionFromMAP: Input is not a valid MAP representation!")    

    return MarginalDistributionFromRAP (D0, D1)

def MarginalDistributionFromMRAP (H):
    """
    Returns the phase type distributed marginal distribution
    of a marked rational arrival process.
    
    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the MRAP
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14
    
    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial vector of the matrix exponentially
        distributed marginal
    A : matrix, shape (M,M)
        The matrix parameter of the matrix exponentially
        distributed marginal    
    """

    if butools.checkInput and not CheckMRAPRepresentation (H):
        raise Exception("MarginalDistributionFromMRAP: Input is not a valid MRAP representation!")    

    Hk = ml.matrix(H[1])
    for i in range (2, len(H)):
        Hk += H[i]
    return (DRPSolve(la.inv(-H[0])*Hk), H[0])

def MarginalDistributionFromMMAP (D):
    """
    Returns the phase type distributed marginal distribution
    of a marked Markovian arrival process.
    
    Parameters
    ----------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the MMAP
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14
    
    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial probability vector of the phase type 
        distributed marginal distribution
    A : matrix, shape (M,M)
        The transient generator of the phase type distributed
        marginal distribution    
    """

    if butools.checkInput and not CheckMMAPRepresentation (D):
        raise Exception("MarginalDistributionFromMMAP: Input is not a valid MMAP representation!")    

    return MarginalDistributionFromMRAP (D)

from butools.ph import MomentsFromPH, MomentsFromME

def MarginalMomentsFromRAP (H0, H1, K=0):
    """
    Returns the moments of the marginal distribution of a 
    rational arrival process.
    
    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational arrival process
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

    if butools.checkInput and not CheckRAPRepresentation (H0, H1):
        raise Exception("MarginalMomentsFromRAP: Input is not a valid RAP representation!")    

    alpha,A = MarginalDistributionFromRAP(H0,H1)
    return MomentsFromME(alpha,A,K)

def MarginalMomentsFromMAP (D0, D1, K=0):
    """
    Returns the moments of the marginal distribution of a 
    Markovian arrival process.
    
    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process
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

    if butools.checkInput and not CheckMAPRepresentation (D0, D1):
        raise Exception("MarginalMomentsFromMAP: Input is not a valid MAP representation!")    

    alpha,A = MarginalDistributionFromMAP(D0,D1)
    return MomentsFromPH(alpha,A,K)
    
def MarginalMomentsFromMRAP (H, K=0):
    """
    Returns the moments of the marginal distribution of a 
    marked rational arrival process.
    
    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the MRAP
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

    if butools.checkInput and not CheckMRAPRepresentation (H):
        raise Exception("MarginalMomentsFromMRAP: Input is not a valid MRAP representation!")    

    alpha,A = MarginalDistributionFromMRAP(H)
    return MomentsFromME(alpha,A,K)

def MarginalMomentsFromMMAP (D, K=0):
    """
    Returns the moments of the marginal distribution of a 
    marked Markovian arrival process.
    
    Parameters
    ----------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the MMAP
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

    if butools.checkInput and not CheckMMAPRepresentation (D):
        raise Exception("MarginalMomentsFromMMAP: Input is not a valid MMAP representation!")    

    alpha,A = MarginalDistributionFromMMAP(D)
    return MomentsFromPH(alpha,A,K)

def LagCorrelationsFromRAP (H0, H1, L=1):
    """
    Returns the lag autocorrelations of a rational arrival
    process.
    
    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational arrival process
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

    if butools.checkInput and not CheckRAPRepresentation (H0, H1):
        raise Exception("LagCorrelationsFromRAP: Input is not a valid RAP representation!")    

    H0i = la.inv(-H0)
    P = H0i*H1
    pi = DRPSolve(P)
    m1, m2 = MomentsFromME(pi, H0, 2)
    pi = pi * H0i * P

    corr = []
    for i in range(L):
        corr.append((np.sum(pi*H0i) - m1*m1) / (m2 - m1*m1))
        pi = pi * P
#    if L>1:
    return np.array(corr)
#    else:
#        return corr[0]

def LagCorrelationsFromMAP (D0, D1, L=1):
    """
    Returns the lag autocorrelations of a Markovian arrival
    process.
    
    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process
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

    if butools.checkInput and not CheckMAPRepresentation (D0, D1):
        raise Exception("LagCorrelationsFromMAP: Input is not a valid MAP representation!")    

    return LagCorrelationsFromRAP (D0, D1, L)

def LagkJointMomentsFromMRAP (H, K=0, L=1):
    """
    Returns the lag-L joint moments of a marked rational 
    arrival process.
    
    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the MRAP to check
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
    Nm : list/cell of matrices of shape(K+1,K+1), length(N)
        The matrices containing the lag-L joint moments,
        starting from moment 0.
    """

    if butools.checkInput and not CheckMRAPRepresentation (H):
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
    pi = DRPSolve(P)
    
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

def LagkJointMomentsFromMMAP (D, K=0, L=1):
    """
    Returns the lag-L joint moments of a marked Markovian 
    arrival process.
    
    Parameters
    ----------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the MMAP to check
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

    if butools.checkInput and not CheckMMAPRepresentation (D):
        raise Exception("LagkJointMomentsFromMMAP: Input is not a valid MMAP representation!")    

    return LagkJointMomentsFromMRAP(D, K, L)

    
def LagkJointMomentsFromRAP (H0, H1, K=0, L=1):
    """
    Returns the lag-L joint moments of a rational arrival
    process.
    
    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational arrival process
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

    if butools.checkInput and not CheckRAPRepresentation (H0, H1):
        raise Exception("LagkJointMomentsFromRAP: Input is not a valid RAP representation!")    

    return LagkJointMomentsFromMRAP((H0,H1), K, L)[0]

def LagkJointMomentsFromMAP (D0, D1, K=0, L=1):
    """
    Returns the lag-L joint moments of a Markovian arrival
    process.
    
    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process
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

    if butools.checkInput and not CheckMAPRepresentation (D0, D1):
        raise Exception("LagkJointMomentsFromMAP: Input is not a valid MAP representation!")    

    return LagkJointMomentsFromRAP(D0, D1, K, L)
    
def RandomMMAP (order, types, mean=1.0, zeroEntries=0, maxTrials=1000, prec=1e-7):
    """
    Returns a random Markovian arrival process with given mean 
    value.
    
    Parameters
    ----------
    order : int
        The size of the MAP
    types : int
        The number of different arrival types
    mean : double, optional
        The mean inter-arrival times of the MMAP
    zeroEntries : int, optional
        The number of zero entries in the D0 and D1 matrices
    maxTrials : int, optional
        The maximum number of trials to find a proper MMAP 
        (that has an irreducible phase process and none of 
        its parameters is all-zero)
    prec : double, optional
        Numerical precision for checking the irreducibility.
        The default value is 1e-14.
    
    Returns
    -------
    D : list/cell of matrices of shape(M,M), length(types+1)
        The D0...Dtypes matrices of the MMAP 
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
                if np.min(np.abs(alpha)) > prec:
                    fullZero = False
                    for Di in D:
                        if np.all(Di==0.0):
                            fullZero = True
                            break
                    if not fullZero:
                        # scale to the mean value
                        m = MarginalMomentsFromMMAP (D, 1)[0]
                        for i in range(types+1):
                            D[i] *= m / mean
                        return D
            trials += 1
    raise Exception("No feasible random MAP/MMAP found with such many zero entries! Try to increase the maxTrials parameter!")

def RandomMAP (order, mean=1.0, zeroEntries=0, maxTrials=1000, prec=1e-7):
    """
    Returns a random Markovian arrival process with given mean 
    value.
    
    Parameters
    ----------
    order : int
        The size of the MAP
    mean : double, optional
        The mean inter-arrival times of the MAP
    zeroEntries : int, optional
        The number of zero entries in the D0 and D1 matrices
    maxTrials : int, optional
        The maximum number of trials to find a proper MAP 
        (that has an irreducible phase process and none of 
        its parameters is all-zero)
    prec : double, optional
        Numerical precision for checking the irreducibility.
        The default value is 1e-14.
    
    Returns
    -------
    D0 : vector, shape (1,M)
        The D0 matrix of the MAP
    D1 : matrix, shape (M,M)
        The D1 matrix of the MAP
    """
    
    return RandomMMAP (order, 1, mean, zeroEntries, maxTrials, prec)
