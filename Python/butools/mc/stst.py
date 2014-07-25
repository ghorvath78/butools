# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 18:29:52 2013

@author: gabor
"""
import numpy as np
import numpy.matlib as ml
import numpy.linalg as la
import butools
from butools.mc import CheckGenerator, CheckProbMatrix


def CRPSolve (Q, prec=1e-14):
    
    if butools.checkInput and np.any(np.sum(Q,1)>prec):
        raise Exception("CRPSolve: The matrix has a rowsum which isn't zero!")
    
    M = np.array(Q)
    M[:,0] = np.ones(M.shape[0])
    m = np.zeros(M.shape[0])
    m[0] = 1.0
    return ml.matrix(la.solve (M.T, m))
    
def CTMCSolve (Q, prec=1e-14):

    if butools.checkInput and not CheckGenerator(Q, False, prec):
        raise Exception("CTMCSolve: The given matrix is not a valid generator. If you are sure you want this use CRPSolve instead of CTMCSolve.");

    return CRPSolve(Q, prec)
    
def DRPSolve (P, prec=1e-14):

    if butools.checkInput and np.any(np.sum(P,1)-1.0>prec):
        raise Exception("DRPSolve: The matrix has a rowsum which isn't 1!");

    return CRPSolve(P-ml.eye(P.shape[0]), prec)
    
def DTMCSolve (P, prec=1e-14):

    if butools.checkInput and not CheckProbMatrix(P, False, prec):
        raise Exception("CTMCSolve: The given matrix is not a valid generator. If you are sure you want this use CRPSolve instead of CTMCSolve.");

    return DRPSolve(P, prec)

