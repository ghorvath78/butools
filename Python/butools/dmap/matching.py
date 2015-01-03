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

    return DMRAPFromMoments (moms, [Nm])
    
def CanonicalFromDMAP2 (D0, D1, prec=1e-14):

    if butools.checkInput:
        if D0.shape[0]!=2:
            raise Exception("CanonicalFromDMAP2: size is not 2!")
        if not CheckDMAPRepresentation(D0, D1, prec):
	        raise Exception("CanonicalFromDMAP2: Input is not a valid DMAP representation!")

    ev = la.eigvals(D0)
    ix = np.argsort(-np.abs(np.real(ev)))
    ev = ev[ix]

    s1=ev[0]
    s2=ev[1]

    if s2>=0:
        G0, G1 = CanonicalFromMAP2 (D0-ml.eye(2),D1,prec)
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

def DMAP2FromMoments (moms, corr1, prec=1e-14):

    Nm = ml.matrix([[1, moms[0]],[moms[0], corr1*(moms[1]-moms[0]**2)+moms[0]**2]])
    
    H0, H1 = DRAPFromMoments (moms, Nm)
   
    oldCheckInput = butools.checkInput
    butools.checkInput = False
    
    D0, D1 = CanonicalFromDMAP2 (H0, H1, prec)
    
    butools.checkInput = oldCheckInput
    
    return (D0,D1)
