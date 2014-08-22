# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 15:28:54 2013

@author: gabor
"""
import numpy as np
import numpy.matlib as ml
import scipy.linalg as la

def TransformToOnes (clovec):
    """
    Returns the similarity transformation matrix that converts 
    the given column vector to a vector of ones. It works even
    if it has zero entries.
    
    Parameters
    ----------
    clovec : column vector, shape(M,1)
        The original closing vector
        
    Returns
    -------
    B : matrix, shape(M,M)
        The matrix by which `B\cdot clovec = \mathbf{1}` holds
    """

    # construct permutation matrix to move at least one non-zero element to the first position
    # to acchieve it, the code below sorts it in a reverse order
    m = clovec.shape[0]    
    ix = np.argsort(-np.array(clovec).flatten())
    P=ml.zeros((m,m))
    for i in range(m):
        P[i,ix[i]] = 1.0
    cp = P*clovec

    # construct transformation matrix B for which B*rp=1 holds
    B = ml.zeros((m,m))
    for i in range(m):
        B[i,0:i+1] = 1.0 / np.sum(cp[0:i+1,0])
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
       
    Ax = ml.matrix(A2)
    Ax[:,0]=Ax[:,0] + ml.ones((N2,1))
    C = ml.zeros((N1,N2))
    C[:,0] = C[:,0] - ml.ones((N1,1))
    return la.solve_sylvester(A1,-Ax,C)
