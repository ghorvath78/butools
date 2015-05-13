# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 15:28:54 2013

@author: gabor
"""
import numpy as np
import numpy.matlib as ml
import scipy.linalg as la
from butools.utils import Linsolve

def SimilarityMatrixForVectors (vecA, vecB):
    """
    Returns the similarity transformation matrix that converts 
    a given column vector to an other column vector. It works 
    even with zero entries.
    
    Parameters
    ----------
    vecA : column vector, shape(M,1)
        The original column vector
    vecB : column vector, shape(M,1)
        The target column vector
        
    Returns
    -------
    B : matrix, shape(M,M)
        The matrix by which `B\cdot vecA = vecB` holds
    """

    # construct permutation matrix to move at least one non-zero element to the first position
    # to acchieve it, the code below sorts it in a reverse order
    m = vecA.shape[0]    
    ix = np.argsort(-np.array(vecA).flatten())
    P=ml.zeros((m,m))
    for i in range(m):
        P[i,ix[i]] = 1.0
    cp = P*vecA

    # construct transformation matrix B for which B*rp=1 holds
    B = ml.zeros((m,m))
    for i in range(m):
        B[i,0:i+1] = vecB.flat[i] / np.sum(cp[0:i+1,0])
    # construct matrix Bp for which Bp*r=1 holds
    return B*P
   
def SimilarityMatrix (A1, A2):
    """
    Returns the matrix that transforms A1 to A2.
    
    Parameters
    ----------
    A1 : matrix, shape (N,N)
        The smaller matrix
    A2 : matrix, shape (M,M)
        The larger matrix (M>=N)
    
    Returns
    -------
    B : matrix, shape (N,M)
        The matrix satisfying `A_1\,B = B\,A_2`
        
    Notes
    -----
    For the existence of a (unique) solution the larger 
    matrix has to inherit the eigenvalues of the smaller one.
    """

    if A1.shape[0]!=A1.shape[1] or A2.shape[0]!=A2.shape[1]:
        raise Exception("SimilarityMatrix: The input matrices must be square!")

    N1 = A1.shape[0]
    N2 = A2.shape[1]

    if N1>N2:
        raise Exception("SimilarityMatrix: The first input matrix must be smaller than the second one!")

    [R1,Q1]=la.schur(A1,'complex')
    [R2,Q2]=la.schur(A2,'complex')
    Q1 = ml.matrix(Q1)
    Q2 = ml.matrix(Q2)
    
    c1 = ml.matrix(np.sum(Q2.H,1))
    c2 = np.sum(Q1.H,1)
    I = ml.eye(N2)
    X = ml.zeros((N1,N2), dtype=complex)
    for k in range(N1-1,-1,-1):
        M = R1[k,k]*I-R2
        if k==N1:
            m = ml.zeros((1,N2))
        else:
            m = -R1[k,k+1:]*X[k+1:,:]
        X[k,:] = Linsolve(np.hstack((M,c1)),np.hstack((m,c2[k])))
    return (Q1*X*ml.matrix(Q2).H ).real
