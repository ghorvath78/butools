# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 17:20:24 2014

@author: gabor
"""

import numpy as np
import numpy.matlib as ml
import scipy.linalg as la
import butools
from butools.mam import *
from butools.mc import CheckGenerator, CheckProbMatrix

def TestMAMPackage ():

    print("---BuTools: MAM package test file---")
    
    print("Enable the verbose messages with the BuToolsVerbose flag")
    butools.verbose = True
    
    print("Enable input parameter checking with the BuToolsCheckInput flag")
    butools.checkInput = True
    
    B = ml.matrix("[0,0;3,4]")
    L = ml.matrix("[-6,5;3,-12]")
    F = ml.matrix("[1,0;2,0]")
    L0 = ml.matrix("[-6,5;6,-8]")
    
    print('----------------------------------------------------------------------------')
    
    print('Input:')
    print('------')
    print("B=\n",B)
    print("L=\n",L)
    print("F=\n",F)
    print("L0=\n",L0)
    
    print('Test:')
    print('-----')
    
    print('R, G, U = QBDFundamentalMatrices (B,L,F,"RGU"):')
    R, G, U = QBDFundamentalMatrices (B,L,F,"RGU")
    print("R=\n",R)
    print("G=\n",G)
    print("U=\n",U)
    
    assert CheckGenerator(U,1)==1, 'QBDFundamentalMatrices: matrix U is not a transient generator!'
    assert np.all(np.abs(la.eigvals(R))<1), 'QBDFundamentalMatrices: the eigenvalues of matrix R are not inside the unit disc!'
    assert CheckProbMatrix(G)==1, 'QBDFundamentalMatrices: matrix G is not a transition prob. matrix!'
    
    print('----------------------------------------------------------------------------')
    
    print('Test:')
    print('-----')
    
    print('pi0, R = QBDSolve (B, L, F, L0):')
    pi0, R = QBDSolve (B, L, F, L0)
    print("pi0=\n",pi0)
    print("R=\n",R)
    
    assert np.sum(pi0)>0 and np.sum(pi0)<=1 and np.all(pi0>=0), 'QBDSolve: wrong pi0 vector!'
    assert np.all(np.abs(la.eigvals(R))<1), 'QBDSolve: the eigenvalues of matrix R are not inside the unit disc!'
    
    
    print('----------------------------------------------------------------------------')
    
    print('Test:')
    print('-----')
    
    print('pi = QBDStationaryDistr (pi0, R, 5):')
    pi = QBDStationaryDistr (pi0, R, 5)
    print(pi)
    
    assert np.sum(pi)>0 and np.sum(pi)<=1 and np.all(pi>=0), 'QBDStationaryDistr: wrong pi vector!'
    
    #==================== M/G/1 type queues ==============================
    B0 = ml.matrix("[0.1, 0.5; 0.3, 0.4]")
    B1 = ml.matrix("[0, 0.1; 0, 0]")
    B2 = ml.matrix("[0.2, 0; 0, 0.2]")
    B3 = ml.matrix("[0, 0.1; 0.1, 0]")
    A0 = ml.matrix("[0.4, 0.2; 0.3, 0.4]")
    A1 = ml.matrix("[0, 0.1; 0, 0]")
    A2 = ml.matrix("[0, 0.2; 0, 0.2]")
    A3 = ml.matrix("[0.1, 0; 0.1, 0]")
    
    B = (B0,B1,B2,B3)
    A = (A0,A1,A2,A3)
    
    G = MG1FundamentalMatrix (A)
    pi = MG1StationaryDistr (A,B,G,300)
    print("G=\n",G)
    print("pi=\n",pi)
    
    assert CheckProbMatrix(G)==1, 'MG1FundamentalMatrix: matrix G is not a transition prob. matrix!'
    assert np.sum(pi)>0 and np.sum(pi)<=1 and np.all(pi>=0), 'MG1StationaryDistr: wrong pi vector!'
    
    #==================== G/M/1 type queues ==============================
    B0 = ml.matrix("[0.7, 0.2; 0.3, 0.6]")
    B1 = ml.matrix("[0.3, 0.4; 0.5, 0.2]")
    B2 = ml.matrix("[0.2, 0.4; 0.1, 0.6]")
    B3 = ml.matrix("[0, 0.1; 0.2, 0]")
    A0 = ml.matrix("[0.1, 0; 0, 0.1]")
    A1 = ml.matrix("[0, 0.2; 0, 0.2]")
    A2 = ml.matrix("[0, 0.1; 0, 0]")
    A3 = ml.matrix("[0.3, 0.2; 0.3, 0.2]")
    A4 = ml.matrix(B3)
    
    B = (B0,B1,B2,B3)
    A = (A0,A1,A2,A3,A4)
    
    R = GM1FundamentalMatrix (A)
    pi = GM1StationaryDistr (B,R,300)
    
    assert np.all(np.abs(la.eigvals(R))<1), 'GM1FundamentalMatrix: the eigenvalues of matrix R are not inside the unit disc!'
    assert np.sum(pi)>0 and np.sum(pi)<=1 and np.all(pi>=0), 'GM1StationaryDistr: wrong pi vector!'
    
    #==================== fluid tests ==============================
    
    Fpp=ml.matrix("[-5 1; 2 -3]")
    Fpm=ml.matrix("[2 1 1; 1 0 0]")
    Fmm=ml.matrix("[-8 4 1; 2 -12 3; 2 0 -2]")
    Fmp=ml.matrix("[3 0; 2 5; 0 0]")
    
    QA = ml.matrix("[-2 2; 5 -5]")
    RA = ml.diag([3, 7])
    NA = QA.shape[0]
    IA = ml.eye(NA)
    QS = ml.matrix("[-4 1 3; 6 -8 2; 3 7 -10]")
    RS = ml.diag([1, 7, 15])
    NS = QS.shape[0]
    IS = ml.eye(NS)
    
    Q = np.kron(QA, IS) + np.kron(IA, QS)
    R = np.kron(RA, IS) - np.kron(IA, RS)
    
    print('----------------------------------------------------------------------------')
    
    print('Input:')
    print('------')
    print("Fpp=",Fpp)
    print("Fpm=",Fpm)
    print("Fmp=",Fmp)
    print("Fmm=",Fmm)
    
    print('Test:')
    print('-----')
    
    print('Psi, K, U = FluidFundamentalMatrices (Fpp, Fpm, Fmp, Fmm, ''PKU''):')
    Psi, K, U = FluidFundamentalMatrices (Fpp, Fpm, Fmp, Fmm, 'PKU')
    print("Psi=",Psi)
    print("K=",K)
    print("U=",U)

    assert CheckGenerator(U,1e-13), 'FluidFundamentalMatrices: matrix U is not a generator!'
    assert np.all(la.eigvals(K)<0), 'FluidFundamentalMatrices: the eigenvalues of matrix K are not negative!'
    assert np.all(Psi>=0) and np.all(Psi<=1) and la.norm(np.sum(Psi,1)-1)<1e-14, 'FluidFundamentalMatrices: matrix Psi is not a transition prob. matrix!'
    
    print('----------------------------------------------------------------------------')
    
    print('Test:')
    print('-----')
    
    print('mass0, ini, K, clo = FluidSolve (Fpp, Fpm, Fmp, Fmm):')
    mass0, ini, K, clo = FluidSolve (Fpp, Fpm, Fmp, Fmm)
    print("mass0=",mass0)
    print("ini=",ini)
    print("K=",K)
    print("clo=",clo)
    
    one = np.sum(mass0) + np.sum(ini*-K.I*clo)
    assert np.abs(one-1)<1e-14, 'FluidSolve: the integral of the fluid distribution in not one!'
    
    print('----------------------------------------------------------------------------')
    
    print('Test:')
    print('-----')
    
    print('[mass0, ini, K, clo] = GeneralFluidSolve (Q, R):')
    mass0, ini, K, clo = GeneralFluidSolve (Q, R)
    print("mass0=",mass0)
    print("ini=",ini)
    print("K=",K)
    print("clo=",clo)
    
    one = np.sum(mass0) + np.sum(ini*-K.I*clo)
    assert np.abs(one-1)<1e-14, 'GeneralFluidSolve: the integral of the fluid distribution in not one!'
    
    print('----------------------------------------------------------------------------')
    
    print('Test:')
    print('-----')
    
    print('y = FluidStationaryDistr (mass0, ini, K, clo, (0:0.1:1)''):')
    y = FluidStationaryDistr (mass0, ini, K, clo, np.linspace(0,30,31))
    print(y)
    
    pi=CTMCSolve(Q)
    assert la.norm(y[-1,:]-pi)<1e-5, 'FluidStationaryDistr: stationary distribution does not converge to the steady state distribution of the phases!'

if __name__ == "__main__":
    TestMAMPackage()