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
from math import sqrt
from butools.ph import MEFromMoments, PH2From3Moments, APHFrom2Moments, APHFrom3Moments
from butools.map import CheckMAPRepresentation, MarginalMomentsFromMAP, LagCorrelationsFromMAP
from butools.moments import MomsFromReducedMoms

def MRAPFromMoments (moms, Nm):

    v, H0 = MEFromMoments (moms)
    H0i = la.inv(-H0)

    Ge = ml.zeros(H0.shape)
    G1 = ml.zeros(H0.shape)

    H0ip = ml.eye(H0.shape[0])
    for i in range(H0.shape[0]):
        Ge[i,:] = v * H0ip
        G1[:,i] = np.sum(H0ip, 1)
        H0ip *= (i+1) * H0i

    Gei = la.inv(Ge)
    G1i = la.inv(G1)
    
    return [-H0*Gei*Nm[i-1]*G1i if i>0 else H0 for i in range(len(Nm)+1)]    

def RAPFromMoments (moms, Nm):

    return MRAPFromMoments (moms, [Nm])
    
def MAP2CorrelationBounds (moms):

    m1, m2, m3 = moms

    h2 = m2/(2.0*m1*m1) - 1
    h3 = m3/(6.0*m1*m1*m1)-m2*m2/(4.0*m1*m1*m1*m1)
    cv2 = float(m2)/m1/m1 - 1.0

    if h2>=0:    
        gub = h2
    else:
        gub = -(h2+math.sqrt(-h3))**2

    if h2<=0 or h3/h2+h2<1:
        glb = -h3 - h2*h2
    else:
        glb = h2 * (h3+h2*h2-h2-math.sqrt((h3+h2*h2-h2)**2+4.0*h2*h2*h2)) / (h3+h2*h2-h2+math.sqrt((h3+h2*h2-h2)**2+4.0*h2*h2*h2))

    if h2>=0:
        return (glb/cv2, gub/cv2)
    else:
        return (gub/cv2, glb/cv2)

def MAP2FromMoments (moms, corr1):

    m1, m2, m3 = moms

    # If we have an exponential distribution, we do not allow correlation
    if abs(m2-2.0*m1*m1) < 1e-14 and abs(corr1) > 1e-14:
        raise Exception ("We do not allow correlation in case of exponentially distributed marginal")

    # Perform PH fitting  
    tau, T =  PH2From3Moments (moms)
    l1 = -T[0,0]
    l2 = -T[1,1]
    p = tau[0,0]
    alpha = l1/l2

    # Check the feasibility of the correlation parameter
    corrl, corru = MAP2CorrelationBounds (moms)
    if corr1 < corrl:
        raise Exception ("The correlation parameter is too small!");
    if corr1 > corru:
        raise Exception ("The correlation parameter is too large!");

    gamma = corr1 * (m2-m1*m1) / (m2/2.0-m1*m1)

    # Perform MAP fitting
    if gamma > 0:
        a = (1.0+alpha*gamma-p*(1.0-gamma)-math.sqrt((1.0+alpha*gamma-p*(1.0-gamma))**2-4.0*alpha*gamma)) / (2.0*alpha)
        b = (1.0+alpha*gamma-p*(1.0-gamma)+math.sqrt((1.0+alpha*gamma-p*(1.0-gamma))**2-4.0*alpha*gamma)) / 2.0
        D0 = ml.matrix ([[-l1, (1.0-a)*l1], [0, -l2]])
        D1 = ml.matrix ([[a*l1, 0], [(1.0-b)*l2, b*l2]])
    elif gamma < 0:
        a = gamma / (alpha*gamma-p*(1.0-gamma))
        b = p*(1.0-gamma)-alpha*gamma
        D0 = ml.matrix ([[-l1, (1.0-a)*l1], [0, -l2]])
        D1 = ml.matrix ([[0, a*l1], [b*l2, (1.0-b)*l2]])
    else:
        D0 = ml.matrix([[-l1, l1], [0, -l2]])
        D1 = ml.matrix([[0, 0], [p*l2, (1.0-p)*l2]])
    return (D0, D1)

def CanonicalFromMAP2 (D0, D1, prec=1e-14):

    if butools.checkInput:
        if D0.shape[0]!=2:
            raise Exception("CanonicalFromMAP2: size is not 2!")
        if not CheckMAPRepresentation(D0, D1, prec):
	        raise Exception("CanonicalFromMAP2: Input is not a valid MAP representation!")

    moms = MarginalMomentsFromMAP (D0, D1, 3, prec)
    corr1 = LagCorrelationsFromMAP (D0, D1, 1, prec)
    return MAP2FromMoments (moms, corr1)


def MAPFromFewMomentsAndCorrelations (moms, corr1=0, r=None):

    m1 = moms[0]
    c2 = moms[1]/moms[0]/moms[0] - 1.0
    if len(moms)>2:
        l3 = moms[2]*moms[0]/moms[1]/moms[1] - 1.0
    else:
        l3 = None
    
    if corr1>=0:
        if r!=None and r<=0 and r>=1:
            raise Exception("Parameter r is out of range")
        if r==None:
            r = 2.0 * corr1 / (1.0 + corr1)
            p1 = (1.0 - (1.0 + corr1) / 2.0) / (1.0 + c2)
            p2 = (1.0 - (1.0 + corr1) / 2.0) * c2 / (1.0 + c2)
        else:
            p1 = (1.0 - corr1 / r) / (1.0 + c2)
            p2 = (1.0 - corr1 / r) * c2 / (1.0 + c2)
        m11 = m1 * (1.0 - sqrt(r))
        m12 = m1 * (1.0 + c2 * sqrt(r))
        if l3==None:
            cv21 = (sqrt(c2)*(1.0+c2)*(1.0+sqrt(r))) / (1.0-sqrt(c2)*(-1.0+sqrt(r))+c2*sqrt(r))
            cv22 = - (c2*(1.0+c2)*(-1.0+r)) / ((1.0+c2*sqrt(r))*(1-sqrt(c2)*(-1.0+sqrt(r))+c2*sqrt(r)))
            m21 = (cv21+1.0) * m11**2
            m22 = (cv22+1.0) * m12**2
            alpha1, A1 = APHFrom2Moments([m11, m21])
            alpha2, A2 = APHFrom2Moments([m12, m22])
        else:
            cv21 = (c2+sqrt(r))/(1.0-sqrt(r))
            cv22 = c2 * (1.0-sqrt(r)) / (1.0+c2*sqrt(r))
            l31 = ((1.0+c2)*l3)/(c2*(1.0-sqrt(r))+sqrt((1.0 +c2*sqrt(r))*c2*(1.0-sqrt(r))))
            l32 = ((1.0+c2)*l3)/((1.0+c2*sqrt(r))+sqrt((1.0+c2*sqrt(r))*c2*(1.0-sqrt(r))))
            m21 = (cv21+1.0)*m11**2
            m22 = (cv22+1.0)*m12**2
            m31 = (l31+1)*m21**2/m11
            m32 = (l32+1)*m22**2/m12
            alpha1, A1 = APHFrom3Moments ([m11, m21, m31])
            alpha2, A2 = APHFrom3Moments ([m12, m22, m32])
    else:
        if c2>=1:
            if r!=None and r<=0 and r>=1/c2:
                raise Exception("Parameter r is out of range")
            if r==None:
                r = -(2.0 * corr1) / (1.0 - c2 * corr1)
                p1 = 0.5 *  (1.0 + (1.0 - c2 * corr1)/2.0)
            else:
                p1 = 0.5 *  (1.0 - corr1/r)                
            p2 = p1
        else:
            if r!=None and r<=0 and r>=1:
                raise Exception("Parameter r is out of range")
            if r==None:
                r = -(2.0 * corr1) / (1.0 - corr1)
                p1 = 0.5 *  (1.0 + (1.0 - corr1)/2.0)
            else:
                p1 = 0.5 *  (1.0 - corr1/r)                
            p2 = p1
        m11 = m1 * (1.0 - sqrt(c2*r))
        m12 = m1 * (1.0 + sqrt(c2*r))
        if l3==None:
            cv21 = c2*(1.0-r)/(1.0-sqrt(c2*r))
            cv22 = c2*(1.0-r)/(1.0+sqrt(c2*r))
            m21 = (cv21+1.0) * m11**2
            m22 = (cv22+1.0) * m12**2
            alpha1, A1 = APHFrom2Moments([m11, m21])
            alpha2, A2 = APHFrom2Moments([m12, m22])
        else:
            cv21 = (c2+sqrt(c2*r))/(1.0-sqrt(c2*r))
            cv22 = (c2-sqrt(c2*r))/(1.0+sqrt(c2*r))
            l31 = 2.0*l3 / (1.0-sqrt(c2*r)+sqrt(1.0-c2*r))
            l32 = 2.0*l3 / (1.0+sqrt(c2*r)+sqrt(1.0-c2*r))
            m21 = (cv21 + 1.0) * m11**2
            m22 = (cv22 + 1.0) * m12**2
            m31 = (l31 + 1.0) * m21**2 / m11
            m32 = (l32 + 1.0) * m22**2 / m12
            alpha1, A1 = APHFrom3Moments ([m11, m21, m31])
            alpha2, A2 = APHFrom3Moments ([m12, m22, m32])
    N1 = A1.shape[0]    
    N2 = A2.shape[0]
    D0 = ml.zeros((N1+N2,N1+N2))
    D0[:N1,:N1] = A1
    D0[N1:,N1:] = A2
    D1 = ml.zeros((N1+N2,N1+N2))
    D1[:N1,:N1] = -np.sum(A1,1)*alpha1*(1.0-p1)
    D1[:N1,N1:] = -np.sum(A1,1)*alpha2*p1
    D1[N1:,:N1] = -np.sum(A2,1)*alpha1*p2
    D1[N1:,N1:] = -np.sum(A2,1)*alpha2*(1.0-p2)
    return D0, D1
    
def RAPFromMomentsAndCorrelations (moms, corr):

    alpha, D0 = MEFromMoments (moms)
    M = alpha.shape[1]
    
    if len(corr) < 2*M-3:
        raise Exception("RAPFromMomentsAndCorrelations: The number of correlations given is less than required the 2n-3!")
    
    rcorr=corr[0:2*M-3]/((moms[1]/2.0-moms[0]**2)/(moms[1]-moms[0]**2))
    rcorr = MomsFromReducedMoms (rcorr)
    tmp,X = MEFromMoments (rcorr.tolist())
    
    N = X.shape[1]

    if N+1 != D0.shape[0]:
        raise Exception("RAPFromMomentsAndCorrelations: Correlation order is different from ME order")

    T1 = ml.zeros((N,N))
    for i in range(N):
        for j in range(i+1):
            T1[i,j] = 1.0

    U1 = ml.zeros((N,N))
    for i in range(N):
        for j in range(i,N):
            U1[i,j] = 1.0 / (N-i)

    T2 = ml.zeros((M,M))
    for i in range(M):
        for j in range(i+1):
            T2[i,j] = 1.0

    U2 = ml.zeros((M,M))
    for i in range(M):
        for j in range(i,M):
            U2[i,j] = 1.0 / (M-i)

    Y = -T1.I*U1*X.I*U1.I*T1
    II = ml.eye(N+1)
    II[1:,1:]=Y
    Y=II
    D1=-D0*U2.I*T2*Y*T2.I*U2
    return (D0, D1)
