# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 16:57:48 2014

@author: gabor
"""

import numpy as np
import numpy.linalg as la
import butools

def CheckGenerator (Q, transient=False, prec=1e-14):

    if Q.shape[0]!=Q.shape[1]:
        if butools.verbose:
            print ("CheckGenerator: Generator is not a square matrix!\n")
        return False

    if np.any(np.diag(Q)>=prec):
        if butools.verbose:
            print ("CheckGenerator: The diagonal of the generator is not negative (precision: {0})!\n".format(prec))
        return False

    N = Q.shape[0]
    odQ = Q<-prec
    for i in range(N):
        odQ[i,i] = 0

    if np.sum(np.any(odQ))>0:
        if butools.verbose:
            print ("CheckGenerator: The generator has negative off-diagonal element (precision: {0})!\n".format(prec))
        return False

    if transient:
        if np.max(np.sum(Q,1))>prec:
            if butools.verbose:
                print ("CheckGenerator: The rowsum of the transient generator is greater than 0 (precision: {0})!\n".format(prec))
            return False

        if np.max(np.real(la.eigvals(Q)))>=prec:
            if butools.verbose:
                print ("CheckGenerator: The transient generator has non-negative eigenvalue (precision: {0})!\n".format(prec))
            return False
    else:
        if np.any(np.abs(np.sum(Q,1))>prec):
            if butools.verbose:
                print ("CheckGenerator: The rowsum of the generator is not 0 (precision: {0})!\n".format(prec))
            return False
    return True

def CheckProbMatrix (P, transient=False, prec=1e-14):

    if P.shape[0]!= P.shape[1]:
        if butools.verbose:
            print ("CheckProbMatrix: the matrix is not a square matrix!\n")
        return False

    if np.min(np.min(P))<-prec:
        if butools.verbose:
            print ("CheckProbMatrix: the matrix has negative element (precision: {0})!\n".format(prec))
        return False

    if transient:
        if np.any(np.sum(P,1)-1.0>prec):
            if butools.verbose:
                print ("CheckProbMatrix: The rowsum of the matrix (transient) is not less or equal than 1 (precision: {0})!\n", prec)
            return False

        if np.max(np.real(la.eigvals(P)))>=1.0-prec:
            if butools.verbose:
                print ("CheckProbMatrix: The real part of the largest eigenvalue of the transient matrix is not less than 1 (precision: {0})!\n".format(prec))
            return False
    else:
        if np.any(np.abs(np.sum(P,1)-1.0)>prec):
            if butools.verbose:
                print ("CheckProbMatrix: The rowsum of the matrix is not 1 (precision: {0})!\n".format(prec))
            return False
    return True

def CheckProbVector (pi, sub=False, prec=1e-14):

    if np.min(pi)<-prec:
        if butools.verbose:
            print ("CheckProbVector: The vector has negative element (precision: {0})!\n".format(prec))
        return False

    if sub:
        if np.sum(pi)>1.0+prec:
            if butools.verbose:
                print ("CheckProbVector: The sum of the substochastic vector is not less than 1 (precision: {0})!".format(prec))
            return False
    else:
        if np.abs(np.sum(pi)-1.0)>prec:
            if butools.verbose:
                print ("CheckProbVector: The sum of the vector is not 1 (precision: {0})!\n".format(prec))
            return False
    
    return True
