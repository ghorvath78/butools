# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 12:12:05 2013

@author: gabor
"""

import numpy as np
import numpy.matlib as ml
import numpy.linalg as la
import butools
import math
from butools.mc import DRPSolve
from butools.map import CanonicalFromMAP2
from butools.moments import JFactorialMomsFromJMoms, FactorialMomsFromMoms
from butools.dph import MGFromMoments
from butools.dmap import CheckDMAPRepresentation

def DMRAPFromMoments (moms, Nm):
    """
    Creates a discrete marked rational arrival process that
    has the same marginal and lag-1 joint moments as given 
    (see [1]_).
    
    Parameters
    ----------
    moms : vector of doubles
        The list of marginal moments. To obtain a discrete 
        marked rational process of order M, 2*M-1 marginal 
        moments are required.
    Nm : list of matrices, shape (M,M)
        The list of lag-1 joint moment matrices. The 
        length of the list determines K, the number of arrival 
        types of the discrete rational process.
    
    Returns
    -------
    H : list of matrices, shape (M,M)
        The H0, H1, ..., HK matrices of the discrete marked
        rational process
    
    References
    ----------
    .. [1] Andras Horvath, Gabor Horvath, Miklos Telek, "A 
           traffic based decomposition of two-class queueing
           networks with priority service," Computer Networks 
           53:(8) pp. 1235-1248. (2009)
    """

    v, H0 = MGFromMoments (moms)
    H0i = la.inv(ml.eye(H0.shape[0])-H0)

    Ge = ml.zeros(H0.shape)
    G1 = ml.zeros(H0.shape)

    H0ip = ml.eye(H0.shape[0])
    for i in range(H0.shape[0]):
        Ge[i,:] = v * H0ip
        G1[:,i] = np.sum(H0ip, 1)
        H0ip *= (i+1) * H0i
        if i>0:
            H0ip *= H0

    Gei = la.inv(Ge)
    G1i = la.inv(G1)
    
    H = [H0]
    for i in range(1,len(Nm)+1):
        Nmi = Nm[i-1]
        row1 = np.array([FactorialMomsFromMoms(Nmi[0,1:].A.flatten())])
        col1 = np.array([FactorialMomsFromMoms(Nmi[1:,0].A.flatten())]).T
        mid = JFactorialMomsFromJMoms(Nmi[1:,1:])
        Nmi = np.bmat([[[[Nmi[0,0]]], row1 ], [col1, mid]])
        H.append((ml.eye(H0.shape[0])-H0)*Gei*Nmi*G1i)
    return H

def DRAPFromMoments (moms, Nm):
    """
    Creates a discrete rational arrival process that has the 
    same marginal and lag-1 joint moments as given (see [1]_).
    
    Parameters
    ----------
    moms : vector of doubles
        The list of marginal moments. To obtain a rational 
        process of order M, 2*M-1 marginal moments are 
        required.
    Nm : matrix, shape (M,M)
        The matrix of lag-1 joint moments. 
    
    Returns
    -------
    H0 : matrix, shape (M,M)
        The H0 matrix of the discrete rational process
    H1 : matrix, shape (M,M)
        The H1 matrix of the discrete rational process
    
    References
    ----------
    .. [1] G Horvath, M Telek, "A minimal representation of 
           Markov arrival processes and a moments matching 
           method," Performance Evaluation 64:(9-12) pp. 
           1153-1168. (2007)       
    """

    return DMRAPFromMoments (moms, [Nm])
    
def CanonicalFromDMAP2 (D0, D1):
    """
    Returns the canonical form of an order-2 discrete Markovian
    arrival process.
    
    Parameters
    ----------
    D0 : matrix, shape (2,2)
        The D0 matrix of the DMAP(2)
    D1 : matrix, shape (2,2)
        The D1 matrix of the DMAP(2)
    prec : double, optional
        Numerical precision to check the input, default 
        value is 1e-14
    
    Returns
    -------
    G0 : matrix, shape (1,2)
        The D0 matrix of the canonical DMAP(2)
    G1 : matrix, shape (2,2)
        The D1 matrix of the canonical DMAP(2)
    """

    if butools.checkInput:
        if D0.shape[0]!=2:
            raise Exception("CanonicalFromDMAP2: size is not 2!")
        if not CheckDMAPRepresentation(D0, D1):
	        raise Exception("CanonicalFromDMAP2: Input is not a valid DMAP representation!")

    ev = la.eigvals(D0)
    ix = np.argsort(-np.abs(np.real(ev)))
    ev = ev[ix]

    s1=ev[0]
    s2=ev[1]

    if s2>=0:
        G0, G1 = CanonicalFromMAP2 (D0-ml.eye(2),D1)
        G0 = G0 + ml.eye(2)
        return (G0, G1)

    #s2 is negative
    av = DRPSolve (la.inv(ml.eye(2)-D0)*D1)

    gamma = la.eigvals(la.inv(ml.eye(2)-D0)*D1)
    ix = np.argsort(-np.abs(np.real(gamma)))
    gamma = gamma[ix]
    gamma = gamma[1]

    w1=1.0/(s1-s2)*(np.sum(D0,1)-s2*ml.ones((2,1)))
    w2=ml.ones((2,1))-w1

    W=np.hstack((w1, w2))
    A=(1.0-s1)*(av*W)
    a1=A[0,0]

    if gamma>=0:
        a=-(1/(2*(-1+s1)*(-1+s1+s2)**2))*(1-4*s1+a1*s1+5*s1**2-a1*s1**2-2*s1**3-2*s2-a1*s2+5*s1*s2-3*s1**2*s2+s2**2+a1*s2**2-s1*s2**2-gamma+3*s1*gamma-a1*s1*gamma-3*s1**2*gamma+a1*s1**2*gamma+s1**3*gamma+s2*gamma+a1*s2*gamma-2*s1*s2*gamma+s1**2*s2*gamma-a1*s2**2*gamma+math.sqrt((-1+s1+s2)**2*((-1+s1**2*(-2+gamma)+gamma+s2*(1+a1-a1*gamma)+s1*(3-a1-s2-2*gamma+a1*gamma))**2-4*(-1+s1)*(-s1**3*(-1+gamma)+a1*(-1+s2)*s2*(-1+gamma)+s1**2*(-2+a1+s2+2*gamma-a1*gamma)+s1*(1-a1-s2-gamma+a1*gamma)))))
        b=1+(a*(-1+s1+s2-s1*s2)*gamma)/((a-1)*(-s1*s2+a*(-1+s1+s2)))

        G0=ml.matrix([[s1+s2, a*(1-s1-s2)], [s1*s2/(a*(s1+s2-1)), 0]])
        G1=ml.matrix([[(1-a)*(1-s1-s2), 0], [b*(1+s1*s2/(a*(1-s1-s2))), (1-b)*(1+s1*s2/(a*(1-s1-s2)))]])
    else:
        #gamma<0
        a=(a1*s1-a1*s1**2+s2-a1*s2-3*s1*s2+2*s1**2*s2-s2**2+a1*s2**2+s1*s2**2+s1*gamma-a1*s1*gamma-2*s1**2*gamma+a1*s1**2*gamma+s1**3*gamma+a1*s2*gamma-a1*s2**2*gamma+math.sqrt(-4*(-1+s1)*s1*s2*(-1+s1+s2)*(a1*(s1-s2)*(-1+gamma)+(-1+s1)*(s2+(-1+s1)*gamma))+(a1*(-s1+s1**2+s2-s2**2)*(-1+gamma)+(-1+s1)*((-1+2*s1)*s2+s2**2+(-1+s1)*s1*gamma))**2))/(2*(-1+s1+s2)*(a1*(s1-s2)*(-1+gamma)+(-1+s1)*(s2+(-1+s1)*gamma)))
        b=-((a*(1-s1)*(1-s2)*gamma)/((a-1)*(-a+a*s1+a*s2-s1*s2)))

        G0=ml.matrix([[s1+s2, a*(1-s1-s2)],[s1*s2/(a*(s1+s2-1)), 0]])
        G1=ml.matrix([[0, (1-a)*(1-s1-s2)],[b*(1-s1*s2/(a*(s1+s2-1))), (1-b)*(1-s1*s2/(a*(s1+s2-1)))]])
    return (G0, G1)

def DMAP2FromMoments (moms, corr1):
    """
    Returns a discrete MAP(2) which has the same 3 marginal
    moments and lag-1 autocorrelation as given.
    
    Parameters
    ----------
    moms : vector, length(3)
        First three marginal moments of the inter-arrival times
    corr1 : double
        The lag-1 autocorrelation of the inter-arrival times
    
    Returns
    -------
    D0 : matrix, shape (2,2)
        The D0 matrix of the discrete MAP(2)
    D1 : matrix, shape (2,2)
        The D1 matrix of the discrete MAP(2)
    
    Notes
    -----
    Raises an exception if the moments are not feasible with
    a DMAP(2). This procedure calls :func:`butools.dmap.DRAPFromMoments`
    followed by :func:`butools.dmap.CanonicalFromDMAP2`.
       
    """

    Nm = ml.matrix([[1, moms[0]],[moms[0], corr1*(moms[1]-moms[0]**2)+moms[0]**2]])
    
    H0, H1 = DRAPFromMoments (moms, Nm)
   
    oldCheckInput = butools.checkInput
    butools.checkInput = False
    
    D0, D1 = CanonicalFromDMAP2 (H0, H1)
    
    butools.checkInput = oldCheckInput
    
    return (D0,D1)
