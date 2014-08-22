# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 17:20:24 2014

@author: gabor
"""

import butools
from butools.reptrans import *
import numpy as np
import numpy.matlib as ml
import scipy.linalg as la
from math import sqrt

def TestRepTransPackage ():

    print('---BuTools: RepTrans package test file---')

    print('Enable the verbose messages with the BuToolsVerbose flag')
    butools.verbose = True
    
    print('Enable input parameter checking with the BuToolsCheckInput flag')
    butools.checkInput = True
   
    eps = np.finfo(float).eps   
   
    print('--SimilarityMatrix------------------------------------------------------------------')
    
    print('Input:')
    print('------')
    
    A1 = ml.matrix([[0.2, 0.8, 0], [1.2, -0.4, 0.1], [-0.2, 0.7, 0.5]])
    T = ml.matrix([[1, 2, -4, 6], [0, 8, -9, 7], [-3, 7, 8, -2]])
    A2 = la.pinv(T)*A1*T
    
    print('Test:')
    print('-----')
    
    print('B=SimilarityMatrix(A1,A2):')
    B=SimilarityMatrix(A1,A2)
    print(B)
    err = la.norm(A1*B-B*A2)
    print(err)
    assert err<sqrt(eps), 'The resulting matrix T does not satisfy A1*T = T*A2!'
    
    print('--TransformToAcyclic----------------------------------------------------------------')
    
    print('Input:')
    print('------')
    
    A = ml.matrix([[-0.8, 0.8, 0], [0.1, -0.3, 0.1], [0.2, 0, -0.5]])
    
    print('Test:')
    print('-----')
    
    print('B=TransformToAcyclic(A):')
    B=TransformToAcyclic(A)
    print(B)
    C=SimilarityMatrix(A,B)
    err = la.norm(A*C-C*B)
    print(err)
    assert err<sqrt(eps), 'The original and the transformed matrix are not similar!'
    
    print('--TransformToMonocyclic-------------------------------------------------------------')
    
    print('Input:')
    print('------')
    
    A = ml.matrix([[-1,0,0],[0,-3,2],[0,-2,-3]])
    
    print('Test:')
    print('-----')
    
    print('B=TransformToMonocyclic(A):')
    B=TransformToMonocyclic(A)
    print(B)
    C=SimilarityMatrix(A,B)
    err = la.norm(A*C-C*B)
    print(err)
    assert err<sqrt(eps), 'The original and the transformed matrix are not similar!'
    
    print('--ExtendToMarkovian-----------------------------------------------------------------')
    
    print('Input:')
    print('------')
    
    print('Original PH:')
    alpha = ml.matrix([[0.2, 0.3, 0.5]])
    print(alpha)
    A = ml.matrix([[-1,0,0],[0,-3,0.6],[0,-0.6,-3]])
    print(A)
    
    print('Transformed to Monocyclic:')
    B=TransformToMonocyclic(A)
    print(B)
    C=SimilarityMatrix(A,B)
    beta = alpha*C
    print(beta)
    
    print('Test:')
    print('-----')
    
    print('[m,M]=ExtendToMarkovian(beta,B):')
    m,M=ExtendToMarkovian(beta,B)
    print(m)
    print(M)
    C=SimilarityMatrix(B,M)
    err = la.norm(B*C-C*M)
    print(err)
    assert err<sqrt(eps), 'The original and the transformed matrix are not similar!'
    assert np.min(m)>-1e-14, 'The initial vector is still not Markovian!'
    
    print('--TransformToOnes-------------------------------------------------------------------')
    
    print('Input:')
    print('------')
    
    print('Original column vector:')
    clovec = ml.matrix([[0.0], [0.3], [-1.5], [0.0]])
    
    print('Test:')
    print('-----')
    
    print('B=TransformToOnes (clovec):')
    B = TransformToOnes (clovec)
    print(B)
    
    assert np.max(np.abs(B*clovec-ml.ones(clovec.shape)))<1e-14, 'The resulting matrix T does not satisfy T*clovec = ones!'
    
    print('--FindMarkovianRepresentation-------------------------------------------------------')    
    print('Tested in the PH and in the MAP package')
    
    print('--MStaircase------------------------------------------------------------------------')   
    print('Tested in the PH and in the MAP package')

if __name__ == "__main__":
    TestRepTransPackage()