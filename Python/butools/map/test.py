# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 17:20:24 2014

@author: gabor
"""

import numpy as np
import numpy.matlib as ml
import butools
from butools.map import *
from butools.moments import CheckMoments
from butools.trace import MarginalMomentsFromTrace
from butools.ph import *

def TestMAPPackage ():

    print("---BuTools: MAP package test file---")
    
    print("Enable the verbose messages with the BuToolsVerbose flag")
    butools.verbose = True
    
    print("Enable input parameter checking with the BuToolsCheckInput flag")
    butools.checkInput = True
    
    print("--MarginalDistributionFromMAP--------------------------------------------------------------------")
    help 
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-0.17, 0, 0, 0.07],[0.01, -0.78, 0.03, 0.08],[0.22, 0.17, -1.1, 0.02],[0.04, 0.12, 0, -0.42]])
    D1=ml.matrix([[0, 0.06, 0, 0.04],[0.04, 0.19, 0.21, 0.22],[0.22, 0.13, 0.15, 0.19],[0.05, 0, 0.17, 0.04]])
    
    print("Test:")
    print("-----")
    
    print("[a,A]=MarginalDistributionFromMAP(D0,D1):")
    [a,A]=MarginalDistributionFromMAP(D0,D1)
    print(a)
    print(A)
    
    assert a.shape[1]==D0.shape[1] and CheckPHRepresentation(a,A), "MarginalDistributionFromMAP returned a wrong PH representation!"
    
    print("--MarginalMomentsFromMAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-0.17, 0, 0, 0.07],[0.01, -0.78, 0.03, 0.08],[0.22, 0.17, -1.1, 0.02],[0.04, 0.12, 0, -0.42]])
    D1=ml.matrix([[0, 0.06, 0, 0.04],[0.04, 0.19, 0.21, 0.22],[0.22, 0.13, 0.15, 0.19],[0.05, 0, 0.17, 0.04]])
    
    print("Test:")
    print("-----")
    
    print("moms=MarginalMomentsFromMAP(D0,D1):")
    moms=MarginalMomentsFromMAP(D0,D1)
    print(moms)
    
    assert len(moms)==2*D0.shape[1]-1 and CheckMoments(moms), "MarginalMomentsFromMAP returned wrong moments!"
    
    print("--MarginalDistributionFromRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[-2, 0, 0],[0, -3, 1],[0, -1, -2]])
    H1=ml.matrix([[1.8, 0.2, 0],[0.2, 1.8, 0],[0.2, 1.8, 1]])
    
    print("Test:")
    print("-----")
    
    print("[a,A]=MarginalDistributionFromRAP(H0,H1):")
    [a,A]=MarginalDistributionFromRAP(H0,H1)
    print(a)
    print(A)
    
    assert a.shape[1]==H0.shape[1] and CheckMERepresentation(a,A), "MarginalDistributionFromRAP returned a wrong ME representation!"
    
    print("--MarginalMomentsFromRAP--------------------------------------------------------------------------")
    help 
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[-2, 0, 0],[0, -3, 1],[0, -1, -2]])
    H1=ml.matrix([[1.8, 0.2, 0],[0.2, 1.8, 0],[0.2, 1.8, 1]])
    
    print("Test:")
    print("-----")
    
    print("moms=MarginalMomentsFromRAP(H0,H1):")
    moms=MarginalMomentsFromRAP(H0,H1)
    print(moms)
    
    assert len(moms)==2*H0.shape[1]-1 and CheckMoments(moms), "MarginalMomentsFromRAP returned wrong moments!"
    
    print("--MarginalDistributionFromMMAP--------------------------------------------------------------------------")
    help 
    
    D0=ml.matrix([[-1.78, 0.29],[0.07, -0.92]])
    D1=ml.matrix([[0.15, 0.49],[0.23, 0.36]])
    D2=ml.matrix([[0.11, 0.2],[0.01, 0]])
    D3=ml.matrix([[0.14, 0.4],[0.11, 0.14]])
    
    print("Test:")
    print("-----")
    
    print("[a,A]=MarginalDistributionFromMMAP((D0,D1,D2,D3)):")
    [a,A]=MarginalDistributionFromMMAP((D0,D1,D2,D3))
    print(a)
    print(A)
    
    assert a.shape[1]==D0.shape[1] and CheckPHRepresentation(a,A), "MarginalDistributionFromMMAP returned a wrong PH representation!"
    
    print("--MarginalMomentsFromMMAP--------------------------------------------------------------------------")
    
    D0=ml.matrix([[-1.78, 0.29],[0.07, -0.92]])
    D1=ml.matrix([[0.15, 0.49],[0.23, 0.36]])
    D2=ml.matrix([[0.11, 0.2],[0.01, 0]])
    D3=ml.matrix([[0.14, 0.4],[0.11, 0.14]])
    
    print("Test:")
    print("-----")
    
    print("moms=MarginalMomentsFromMMAP((D0,D1,D2,D3)):")
    moms=MarginalMomentsFromMMAP((D0,D1,D2,D3))
    print(moms)
    
    assert len(moms)==2*D0.shape[1]-1 and CheckMoments(moms), "MarginalMomentsFromMMAP returned wrong moments!"
    
    print("--MarginalDistributionFromMRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    x=0.18
    H0=ml.matrix([[-5, 0.1+x, 0.9, 1],[1, -8, 0.9, 0.1],[0.9, 0.1, -4, 1],[1, 2, 3, -9]])
    H1=ml.matrix([[0.1-x, 0.7, 0.1, 0.1],[0.1, 1, 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
    H2=ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1, 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1, 0.1, 0.8]])
    
    print("Test:")
    print("-----")
    
    print("[a,A]=MarginalDistributionFromMRAP((H0,H1,H2)):")
    [a,A]=MarginalDistributionFromMRAP((H0,H1,H2))
    print(a)
    print(A)
    
    assert a.shape[1]==H0.shape[1] and CheckMERepresentation(a,A), "MarginalDistributionFromMRAP returned a wrong ME representation!"
    
    print("--MarginalMomentsFromMRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    x=0.18
    H0=ml.matrix([[-5, 0.1+x, 0.9, 1],[1, -8, 0.9, 0.1],[0.9, 0.1, -4, 1],[1, 2, 3, -9]])
    H1=ml.matrix([[0.1-x, 0.7, 0.1, 0.1],[0.1, 1, 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
    H2=ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1, 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1, 0.1, 0.8]])
    
    print("Test:")
    print("-----")
    
    print("moms=MarginalMomentsFromMRAP((H0,H1,H2)):")
    moms=MarginalMomentsFromMRAP((H0,H1,H2))
    print(moms)
    
    assert len(moms)==2*H0.shape[1]-1 and CheckMoments(moms), "MarginalMomentsFromMRAP returned wrong moments!"
    
    print("--LagCorrelationsFromMAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-5, 0, 1, 1],[1, -8, 1, 0],[1, 0, -4, 1],[1, 2, 3, -9]])
    D1=ml.matrix([[0, 1, 0, 2],[2, 1, 3, 0],[0, 0, 1, 1],[1, 1, 0, 1]])
    
    print("Test:")
    print("-----")
    
    print("LagCorrelationsFromMAP(D0,D1,3):")
    corr = LagCorrelationsFromMAP(D0,D1,3)
    print(corr)
    
    assert len(corr)==3 and np.all(np.array(corr)<1) and np.all(np.array(corr)>-1), "LagCorrelationsFromMAP returned wrong autocorrelation coefficients!"
    
    print("--LagCorrelationsFromRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[-2, 0, 0],[0, -3, 1],[0, -1, -2]])
    H1=ml.matrix([[1.8, 0.2, 0],[0.2, 1.8, 0],[0.2, 1.8, 1]])
    
    print("Test:")
    print("-----")
    
    print("LagCorrelationsFromRAP(H0,H1,3):")
    corr = LagCorrelationsFromRAP(H0,H1,3)
    print(corr)
    
    assert len(corr)==3 and np.all(np.array(corr)<1) and np.all(np.array(corr)>-1), "LagCorrelationsFromRAP returned wrong autocorrelation coefficients!"
    
    print("--LagkJointMomentsFromMAP--------------------------------------------------------------------------")
    help 
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-5, 0, 1, 1],[1, -8, 1, 0],[1, 0, -4, 1],[1, 2, 3, -9]])
    D1=ml.matrix([[0, 1, 0, 2],[2, 1, 3, 0],[0, 0, 1, 1],[1, 1, 0, 1]])
    
    print("Test:")
    print("-----")
    
    print("LagkJointMomentsFromMAP(D0,D1,2,1):")
    Nm=LagkJointMomentsFromMAP(D0,D1,4,1)
    print(Nm)
    
    moms=MarginalMomentsFromMAP(D0,D1,4)
    print(moms)
    
    print("Correlation from joint moments:")
    cjm=np.empty(3)
    for i in range(3):
        Nx=LagkJointMomentsFromMAP(D0,D1,1,i+1)
        print(Nx)
        print(Nx[1,1])
        cjm[i] = (Nx[1,1]-moms[0]**2) / (moms[1]-moms[0]**2)
    print(cjm)
    corr = LagCorrelationsFromMAP(D0,D1,3)
    print(corr)
    
    assert np.all(np.all(Nm>0)) and la.norm(moms-Nm[0,1:])<1e-14 and la.norm(moms-Nm[1:,0].A.flatten())<1e-14 and la.norm(corr-cjm)<1e-14, "Joint moment matrix is invalid!"
    

    print("--LagkJointMomentsFromRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[-2, 0, 0],[0, -3, 1],[0, -1, -2]])
    H1=ml.matrix([[1.8, 0.2, 0],[0.2, 1.8, 0],[0.2, 1.8, 1]])
    
    print("Test:")
    print("-----")
    
    print("LagkJointMomentsFromRAP(H0,H1,2,1):")
    Nm=LagkJointMomentsFromRAP(H0,H1,4,1)
    print(len(Nm))
    
    
    moms=MarginalMomentsFromRAP(H0,H1,4)
    
    print("Correlation from joint moments:")
    cjm=np.empty(3)
    for i in range(3):
        Nx=LagkJointMomentsFromRAP(H0,H1,1,i+1)
        print(type(Nx))
        cjm[i] = (Nx[1,1]-moms[0]**2) / (moms[1]-moms[0]**2)
    print(cjm)
    corr = LagCorrelationsFromRAP(H0,H1,3)
    print(corr)
    
    assert np.all(np.all(Nm>0)) and la.norm(moms-Nm[0,1:])<1e-14 and la.norm(moms-Nm[1:,0].A.flatten())<1e-14 and la.norm(corr-cjm)<1e-14, "Joint moment matrix is invalid!"
    
    print("--LagkJointMomentsFromMMAP--------------------------------------------------------------------------")
    help 
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-1.78, 0.29],[0.07, -0.92]])
    D1=ml.matrix([[0.15, 0.49],[0.23, 0.36]])
    D2=ml.matrix([[0.11, 0.2],[0.01, 0]])
    D3=ml.matrix([[0.14, 0.4],[0.11, 0.14]])
    
    print("Test:")
    print("-----")
    
    print("LagkJointMomentsFromMMAP((D0,D1,D2,D3),3,1)")
    Nm=LagkJointMomentsFromMMAP((D0,D1,D2,D3),3,1)
    
    print("Moments of arrival type 1, Nm[0]:")
    print(Nm[0])
    print("Moments of arrival type 2, Nm[1]:")
    print(Nm[1])
    print("Moments of arrival type 2, Nm[2]:")
    print(Nm[2])
    
    assert len(Nm)==3 and la.norm(SumMatrixList(Nm)-LagkJointMomentsFromMAP(D0,D1+D2+D3,3,1))<1e-14,"Joint moment matrix is invalid!"
    
    print("--LagkJointMomentsFromMRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    x=0.18
    H0=ml.matrix([[-5, 0.1+x, 0.9, 1],[1, -8, 0.9, 0.1],[0.9, 0.1, -4, 1],[1, 2, 3, -9]])
    H1=ml.matrix([[0.1-x, 0.7, 0.1, 0.1],[0.1, 1, 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
    H2=ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1, 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1, 0.1, 0.8]])

    print("Test:")
    print("-----")
    
    print("LagkJointMomentsFromMRAP((H0,H1,H2),3,2)")
    Nm=LagkJointMomentsFromMRAP((H0,H1,H2),3,2)
    
    print("Moments of arrival type 1, Nm[0]:")
    print(Nm[0])
    print("Moments of arrival type 2, Nm[1]:")
    print(Nm[1])
    
    assert len(Nm)==2 and la.norm(SumMatrixList(Nm)-LagkJointMomentsFromRAP(H0,H1+H2,3,2))<1e-14,"Joint moment matrix is invalid!"
    
    print("--RandomMAP--------------------------------------------------------------------------")
    
    print("Test:")
    print("-----")
    
    print("[D0,D1]=RandomMAP(4,1.62,10):")
    [D0,D1]=RandomMAP(4,1.62,10)
    print(D0)
    print(D1)
    
    print("Check the mean of the obtained MAP:")
    m = MarginalMomentsFromMAP(D0,D1,1)[0]
    print(m)
    
    assert CheckMAPRepresentation(D0,D1), "RandomMAP failed to return a valid MAP representation!"
    assert np.max(np.abs(m-1.62))<1e-14, "RandomMAP failed to match the given mean value!"
    
    print("--RandomMMAP--------------------------------------------------------------------------")
    
    print("Test:")
    print("-----")
    
    print("D=RandomMMAP(4,3,1.62,10):")
    D=RandomMMAP(4,3,1.62,10)
    print(D[0])
    print(D[1])
    print(D[2])
    print(D[3])
    
    print("Check the mean of the obtained MMAP:")
    m = MarginalMomentsFromMMAP(D,1)[0]
    print(m)
    
    assert CheckMMAPRepresentation(D), "RandomMMAP failed to return a valid MMAP representation!"
    assert np.max(np.abs(m-1.62))<1e-14, "RandomMMAP failed to match the given mean value!"
    
    print("-CheckMAPRepresentation---------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-1, 0, 1],[0, -2, 0],[1, 0, -3]])
    D1=ml.matrix([[-1, 0, 1, 0],[ 0, -2, 0, 1],[ 1, 0, -3, 0],[1, 2, 2, 1]])
    
    print("Test:")
    print("-----")
    
    print("CheckMAPRepresentation(D0,D1):")
    flag=CheckMAPRepresentation(D0,D1)
    print(flag)
    
    assert flag==False, "CheckMAPRepresentation failed to detect the incompatible shapes of D0 and D1!"
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-1, 0, 1],[0, -2, 0],[1, 0, -3]])
    D1=ml.matrix([[1, 0, 1],[0, 2, 0],[1, 0, 3]])
    
    print("Test:")
    print("-----")
    
    print("CheckMAPRepresentation(D0,D1):")
    flag=CheckMAPRepresentation(D0,D1)
    print(flag)
    
    assert flag==False, "CheckMAPRepresentation failed to detect invalid rowsums!"
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-3, 0, 1],[0, -2, 0],[1, 0, -5]])
    D1=ml.matrix([[1, 0, 1],[0, 2, 0],[1, 0, 3]])
    
    print("Test:")
    print("-----")
    
    print("CheckMAPRepresentation(D0,D1):")
    flag=CheckMAPRepresentation(D0,D1)
    print(flag)
    
    assert flag==True, "CheckMAPRepresentation failed to recognize a valid MAP representation!"
    
    print("--CheckRAPRepresentation--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[-1, 0, 1],[0, -2, 0],[1, 0 ,-3],[1, 2, 2]])
    H1=ml.matrix([[-1, 0, 1],[0, -2, 0],[1, 0, -3],[1, 2, 2]])
    
    print("Test:")
    print("-----")
    
    print("CheckRAPRepresentation(H0,H1):")
    flag=CheckRAPRepresentation(H0,H1)
    print(flag)
    
    assert flag==False, "CheckRAPRepresentation failed to detect the incompatible shapes of D0 and D1!"
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[-1, 0, 2],[0, 2, 0],[1, 0 ,-3]])
    H1=ml.matrix([[-1, 0, 1],[0, -2, 0],[1, 0, -3]])
    
    print("Test:")
    print("-----")
    
    print("CheckRAPRepresentation(H0,H1):")
    flag=CheckRAPRepresentation(H0,H1)
    print(flag)
    
    assert flag==False, "CheckRAPRepresentation failed to detect invalid rowsums!"
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[-1, 0, 0],[0, -2, 2],[0, 3, -3]])
    H1=ml.matrix([[0, 0, 1],[0, -1, 1],[1, 0, -1]])
    
    print("Test:")
    print("-----")
    
    print("CheckRAPRepresentation(H0,H1):")
    flag=CheckRAPRepresentation(H0,H1)
    print(flag)
    
    assert flag==False, "CheckRAPRepresentation failed to detect invalid eigenvalues!"
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[-2, 0, 0],[0, -1, 1],[0, -1, -1]])
    H1=ml.matrix([[1, 0, 1],[0, 1, -1],[1, 0, 1]])

    print("Test:")
    print("-----")
    
    print("CheckRAPRepresentation(H0,H1):")
    flag=CheckRAPRepresentation(H0,H1)
    print(flag)
    
    assert flag==False, "CheckRAPRepresentation failed to detect invalid eigenvalues!"
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[-1, 0, 0],[0, -2, 1],[0, -1, -2]])
    H1=ml.matrix([[1, 0, 0],[0, 1, 0],[1, 1, 1]])
    
    print("Test:")
    print("-----")
    
    print("CheckRAPRepresentation(H0,H1):")
    flag=CheckRAPRepresentation(H0,H1)
    print(flag)
    
    assert flag==True, "CheckRAPRepresentation failed to recognize a valid RAP representation!"
    
    print("-CheckMMAPRepresentation---------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-1.05, 0.03, 0.07],[ 0.19, -1.63, 0.06],[ 0, 0.2, -1.03]])
    D1=ml.matrix([[0.16, 0.11, 0],[ 0.1, 0.16, 0],[ 0.27, 0, 0.19]])
    D2=ml.matrix([[0.01, 0.09, 0.13],[ 0.26, 0.21, 0.05],[ 0, 0.16, 0.07]])
    D3=ml.matrix([[0.19, 0.06, 0.2],[ 0.17, 0.16, 0.27],[ 0, 0, 0.14]])
    
    print("Test:")
    print("-----")
    
    print("CheckMMAPRepresentation((D0,D1,D2,D3)):")
    flag=CheckMMAPRepresentation((D0,D1,D2,D3))
    print(flag)
    
    assert flag==True, "CheckMMAPRepresentation failed to recognize a valid MMAP representation!"
    
    print("-CheckMRAPRepresentation---------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    x=0.18
    H0=ml.matrix([[-5, 0.1+x, 0.9, 1],[1, -8, 0.9, 0.1],[0.9, 0.1, -4, 1],[1, 2, 3, -9]])
    H1=ml.matrix([[0.1-x, 0.7, 0.1, 0.1],[0.1, 1, 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
    H2=ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1, 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1, 0.1, 0.8]])
    
    print("Test:")
    print("-----")
    
    print("CheckMRAPRepresentation((H0,H1,H2)):")
    flag=CheckMRAPRepresentation((H0,H1,H2))
    print(flag)
    
    assert flag==True, "CheckMRAPRepresentation failed to recognize a valid MRAP representation!"
    
    print("--RAPFromMoments--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    G0=ml.matrix([[-6.2, 2, 0],[ 2, -9, 1],[ 1, 0, -3]])
    G1=ml.matrix([[2.2, -2, 4],[ 2, 2, 2],[ 1, 0, 1]])
    moms=MarginalMomentsFromRAP(G0,G1,5)
    Nm=LagkJointMomentsFromRAP(G0,G1,2,1)
    print(moms)
    print(Nm)
    
    print("Test:")
    print("-----")
    
    print("[H0,H1]=RAPFromMoments(moms,Nm):")
    [H0,H1]=RAPFromMoments(moms,Nm)
    print(H0)
    print(H1)
    
    rmoms=MarginalMomentsFromRAP(H0,H1,5)
    rNm=LagkJointMomentsFromRAP(H0,H1,2,1)
    print(rmoms)
    print(rNm)
    
    assert la.norm(np.array(moms)-np.array(rmoms))<1e-12 and la.norm(Nm-rNm)<1e-12, "The moments and joint moments returned by RAPFromMoments are not the same as given!"
    
    print("Input:")
    print("------")
    
    G0=ml.matrix([[-5, 0, 1, 1],[1, -8, 1, 0],[1, 0, -4, 1],[1, 2, 3, -9]])
    G1=ml.matrix([[0, 1, 0, 2],[2, 1, 3, 0],[0, 0, 1, 1],[1, 1, 0, 1]])

    moms=MarginalMomentsFromRAP(G0,G1,7)
    Nm=LagkJointMomentsFromRAP(G0,G1,3,1)
    print(moms)
    print(Nm)
    
    print("Test:")
    print("-----")
    
    print("[H0,H1]=RAPFromMoments(moms,Nm):")
    [H0,H1]=RAPFromMoments(moms,Nm)
    print(H0)
    print(H1)
    
    butools.checkPrecision = 1e-8
    rmoms=MarginalMomentsFromRAP(H0,H1,7)
    rNm=LagkJointMomentsFromRAP(H0,H1,3,1)
    print(rmoms)
    print(rNm)
    
    assert CheckRAPRepresentation(H0,H1,1e-8), "RAPFromMoments returned an invalid RAP representation!"
    assert la.norm(np.array(moms)-np.array(rmoms))<1e-8 and la.norm(Nm-rNm)<1e-8, "The moments and joint moments returned by RAPFromMoments are not the same as given!"
    
    print("--MRAPFromMoments--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    G0=ml.matrix([[-1.05, 0.03, 0.07],[0.19, -1.63, 0.06],[0, 0.2, -1.03]])
    G1=ml.matrix([[0.16, 0.11, 0],[0.1, 0.16, 0],[0.27, 0, 0.19]])
    G2=ml.matrix([[0.01, 0.09, 0.13],[0.26, 0.21, 0.05],[0, 0.16, 0.07]])
    G3=ml.matrix([[0.19, 0.06, 0.2],[0.17, 0.16, 0.27],[0, 0, 0.14]])
    G=(G0,G1,G2,G3)
    moms=MarginalMomentsFromMRAP(G,5)
    print(moms)
    Nm=LagkJointMomentsFromMRAP(G,2,1)
    Nm1, Nm2, Nm3 = Nm
    print(Nm1)
    print(Nm2)
    print(Nm3)
    
    print("Test:")
    print("-----")
    
    print("H=MRAPFromMoments(moms,Nm):")
    H=MRAPFromMoments(moms,Nm)
    print(H[0])
    print(H[1])
    print(H[2])
    print(H[3])
    
    butools.checkPrecision = 1e-11
    rmoms=MarginalMomentsFromMRAP(H,5)
    print(rmoms)
    rNm1, rNm2, rNm3 = LagkJointMomentsFromMRAP(H,2,1)
    print(rNm1)    
    print(rNm2)    
    print(rNm3)    
    
    assert la.norm(np.array(moms)-np.array(rmoms))<1e-9 and la.norm(Nm1-rNm1)<1e-10 and la.norm(Nm2-rNm2)<1e-10 and la.norm(Nm3-rNm3)<1e-10, "The moments and joint moments returned by MRAPFromMoments are not the same as given!"
    
    print("--RAPFromMomentsAndCorrelations--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    H1=ml.matrix([[2.2, 0, 2],[0, 4, 2],[0, 1, 1]])
    mom=MarginalMomentsFromRAP(H0,H1)
    corr=LagCorrelationsFromRAP(H0,H1,3)
    print(mom)
    print(corr)
    
    print("Test:")
    print("-----")
    
    print("[G0,G1]=RAPFromMomentsAndCorrelations(mom,corr)")
    [G0,G1]=RAPFromMomentsAndCorrelations(mom,corr)
    print(G0)
    print(G1)
    
    rmom=MarginalMomentsFromRAP(G0,G1)
    rcorr=LagCorrelationsFromRAP(G0,G1,3)
    print(rmom)
    print(rcorr)
    
    assert CheckRAPRepresentation(G0,G1), "RAPFromMomentsAndCorrelations returned an invalid RAP representation!"
    assert la.norm(np.array(rmom)-np.array(mom))+la.norm(np.array(rcorr)-np.array(corr))<1e-12, "The result of RAPFromMomentsAndCorrelations has different moments or correlations than given!"
    
    print("--MAP2FromMoments--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-14, 1],[1, -25]])
    D1=ml.matrix([[6, 7],[3, 21]])
    moms=MarginalMomentsFromMAP(D0, D1, 3)
    corr=LagCorrelationsFromMAP(D0, D1, 1)
    print(moms)
    print(corr)
    
    print("Test:")
    print("-----")
    
    print("[D0,D1]=MAP2FromMoments(moms,corr):")
    [D0,D1]=MAP2FromMoments(moms,corr)
    print(D0)
    print(D1)
    
    rmoms=MarginalMomentsFromMAP(D0, D1, 3)
    rcorr=LagCorrelationsFromMAP(D0, D1, 1)
    print(rmoms)
    print(rcorr)    
    
    assert CheckMAPRepresentation(D0,D1), "MAP2FromMoments returned an invalid MAP representation!"
    assert la.norm(np.array(moms)-np.array(rmoms))<1e-12 and abs(corr-rcorr)<1e-12, "The moments and the correlation returned by MAP2FromMoments are not the same as given!"
    
    print("--MAP2CorrelationBounds--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-14, 1],[1, -25]])
    D1=ml.matrix([[6, 7],[3, 21]])
    moms=MarginalMomentsFromMAP(D0,D1,3)
    print(moms)
    
    print("Test:")
    print("-----")
    
    print("[lb,ub]=MAP2CorrelationBounds(moms):")
    [lb,ub]=MAP2CorrelationBounds(moms)
    print(lb)
    print(ub)
    
    assert lb<=0 and lb>=-1 and ub>=0 and ub<=1, "Correlation bounds given by MAP2CorrelationBounds are not correct"
    
    print("--MAPFromFewMomentsAndCorrelations--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    moms = [1.1, 6.05]
    corr1 = -0.17
    print(moms)
    print(corr1)
    
    print("Test:")
    print("-----")
    
    print("[D0,D1]=MAPFromFewMomentsAndCorrelations(moms, corr1):")
    [D0,D1]=MAPFromFewMomentsAndCorrelations(moms, corr1)
    print(D0)
    print(D1)
    
    rmoms = MarginalMomentsFromMAP(D0,D1,2)
    rcorr1 = LagCorrelationsFromMAP(D0,D1,1)
    print(rmoms)
    print(rcorr1)
    
    assert CheckMAPRepresentation(D0,D1), "MAPFromFewMomentsAndCorrelations returned with a non-Markovian representation!"
    assert la.norm(np.array(rmoms)-np.array(moms))<1e-12 and abs(rcorr1-corr1)<1e-12, "MAPFromFewMomentsAndCorrelations failed to match the marginal moments or the lag-1 autocorrelation!"
    
    print("Input:")
    print("------")
    
    moms = [1.2, 4.32, 20]
    corr1 = 0.4
    print(moms)
    print(corr1)
    
    print("Test:")
    print("-----")
    
    print("[D0,D1]=MAPFromFewMomentsAndCorrelations(moms, corr1):")
    [D0,D1]=MAPFromFewMomentsAndCorrelations(moms, corr1)
    print(D0)
    print(D1)
    
    butools.checkPrecision = 1e-12
    rmoms = MarginalMomentsFromMAP(D0,D1,3)
    rcorr1 = LagCorrelationsFromMAP(D0,D1,1)
    print(rmoms)
    print(rcorr1)
    
    assert CheckMAPRepresentation(D0,D1,1e-13), "MAPFromFewMomentsAndCorrelations returned with a non-Markovian representation!"
    assert la.norm(np.array(rmoms)-np.array(moms))<1e-12 and abs(rcorr1-corr1)<1e-12, "MAPFromFewMomentsAndCorrelations failed to match the marginal moments or the lag-1 autocorrelation!"
    
    print("Input:")
    print("------")
    
    moms = [1.2, 4.32, 20]
    corr1 = 0.4
    print(moms)
    print(corr1)
    
    print("Test:")
    print("-----")
    
    print("[D0,D1]=MAPFromFewMomentsAndCorrelations(moms, corr1):")
    [D0,D1]=MAPFromFewMomentsAndCorrelations(moms, corr1)
    print(D0)
    print(D1)
    
    rmoms = MarginalMomentsFromMAP(D0,D1,3)
    rcorr1 = LagCorrelationsFromMAP(D0,D1,1)
    print(rmoms)
    print(rcorr1)
    
    assert CheckMAPRepresentation(D0,D1,1e-13), "MAPFromFewMomentsAndCorrelations returned with a non-Markovian representation!"
    assert la.norm(np.array(rmoms)-np.array(moms))<1e-12 and abs(rcorr1-corr1)<1e-12, "MAPFromFewMomentsAndCorrelations failed to match the marginal moments or the lag-1 autocorrelation!"
    
    print("--CanonicalFromMAP2--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-14, 1],[1, -25]])
    D1=ml.matrix([[6, 7],[3, 21]])
    
    print("Test:")
    print("-----")
    
    print("[H0,H1]=CanonicalFromMAP2(D0,D1):")
    [H0,H1]=CanonicalFromMAP2(D0,D1)
    print(H0)
    print(H1)
    
    C=SimilarityMatrix(H0,D0)
    err = la.norm(H0*C-C*D0) + la.norm(H1*C-C*D1)
    print(err)
    
    assert CheckMAPRepresentation(H0,H1), "The result of CanonicalFromMAP2 is not a valid MAP representation!"
    assert err<1e-12, "The MAP returned by CanonicalFromMAP2 is not similar to the input!"
    
    print("--MAPFromRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-2, 2],[2, -9]])
    D1=ml.matrix([[-2, 2],[3, 4]])
    
    print("Test:")
    print("-----")
    
    print("[H0,H1]=MAPFromRAP(D0,D1):")
    [H0,H1]=MAPFromRAP(D0,D1)
    print(H0)
    print(H1)
    
    err = la.norm(LagkJointMomentsFromRAP(D0,D1,3,1)-LagkJointMomentsFromRAP(H0,H1,3,1))
    print(err)
    assert err<1e-12, "The RAP returned by MAPFromRAP is not similar to the input!"
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-2.4, 2],[2, -9]])
    D1=ml.matrix([[-1.6, 2],[3, 4]])
    
    print("Test:")
    print("-----")
    
    print("[H0,H1]=MAPFromRAP(D0,D1):")
    [H0,H1]=MAPFromRAP(D0,D1)
    print(H0)
    print(H1)
    
    err = la.norm(LagkJointMomentsFromRAP(D0,D1,3,1)-LagkJointMomentsFromRAP(H0,H1,3,1))
    print(err)
    assert err<1e-12, "The MAP returned by MAPFromRAP is not similar to the input!"
    assert CheckMAPRepresentation(H0,H1), "The result of MAPFromRAP is not a MAP, as it should be!"
    
    print("--MMAPFromMRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    

    x=0.18
    H0=ml.matrix([[-5, 0.1+x, 0.9, 1],[1, -8, 0.9, 0.1],[0.9, 0.1, -4, 1],[1, 2, 3, -9]])
    H1=ml.matrix([[0.1-x, 0.7, 0.1, 0.1],[0.1, 1, 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
    H2=ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1, 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1, 0.1, 0.8]])
    H=(H0,H1,H2)
    
    moms=MarginalMomentsFromMRAP(H)
    jmom=LagkJointMomentsFromMRAP(H,3,1)
    print(moms)   
    print(jmom)
    
    print("Test:")
    print("-----")
    
    print("G=MMAPFromMRAP({H0,H1,H2}):")
    G=MMAPFromMRAP(H)
    print(G[0])
    print(G[1])
    print(G[2])
    
    rmoms=MarginalMomentsFromMMAP(G)
    rjmom=LagkJointMomentsFromMMAP(G,3,1)
    print(rmoms)   
    print(rjmom)
    
    err = la.norm(rjmom[0]-jmom[0]) + la.norm(rjmom[1]-jmom[1])
    print(err)
    assert err<1e-12, "The MMAP returned by MMAPFromMRAP is not similar to the input!"
    assert CheckMMAPRepresentation(G), "The result of MMAPFromMRAP is not a MMAP, as it should be!"
    
    print("-MinimalRepFromRAP---------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-5, 1, 0],[3, -3, 0],[1, 1, -5]])
    D1=ml.matrix([[0, 0, 4],[0, 0, 0],[1, 1, 1]])
    
    print("Test:")
    print("-----")
    
    print("[H0,H1]=MinimalRepFromRAP(D0,D1,""cont"")")
    [H0,H1]=MinimalRepFromRAP(D0,D1,"cont")
    print(H0)
    print(H1)
    
    C = SimilarityMatrix(H0,D0)
    err = la.norm(H0*C-C*D0) + la.norm(H1*C-C*D1)
    print(err)
    
    assert CheckRAPRepresentation(H0,H1), "MinimalRepFromRAP did not return a valid RAP representation!"
    assert H0.shape[0]==3 and err<1e-12, "MinimalRepFromRAP returned a RAP which is non-similar to the input or has an unexpected size!"
    
    print("[H0,H1]=MinimalRepFromRAP(D0,D1,""obs"")")
    [H0,H1]=MinimalRepFromRAP(D0,D1,"obs")
    print(H0)
    print(H1)
    
    C = SimilarityMatrix(H0,D0)
    err = la.norm(H0*C-C*D0) + la.norm(H1*C-C*D1)
    print(err)
    
    assert CheckRAPRepresentation(H0,H1), "MinimalRepFromRAP did not return a valid RAP representation!"
    assert H0.shape[0]==2 and err<1e-12, "MinimalRepFromRAP returned a RAP which is non-similar to the input or has an unexpected size!"
    
    print("[H0,H1]=MinimalRepFromRAP(D0,D1,""obscont"")")
    [H0,H1]=MinimalRepFromRAP(D0,D1,"obscont")
    print(H0)
    print(H1)
    
    C = SimilarityMatrix(H0,D0)
    err = la.norm(H0*C-C*D0) + la.norm(H1*C-C*D1)
    print(err)
    
    assert CheckRAPRepresentation(H0,H1), "MinimalRepFromRAP did not return a valid RAP representation!"
    assert H0.shape[0]==2 and err<1e-12, "MinimalRepFromRAP returned a RAP which is non-similar to the input or has an unexpected size!"
    
    print("--MinimalRepFromMRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-5, 1, 0],[3, -3, 0],[1, 1, -5]])
    D1=0.2*ml.matrix([[0, 0, 4],[0, 0, 0],[1, 1, 1]])
    D2=0.8*ml.matrix([[0, 0, 4],[0, 0, 0],[1, 1, 1]])
    D=(D0,D1,D2)
    
    print("Test:")
    print("-----")
    
    print("H=MinimalRepFromMRAP(D,""cont"")")
    H=MinimalRepFromMRAP(D,"cont")
    print(H[0])
    print(H[1])
    print(H[2])
    
    C = SimilarityMatrix(H[0],D[0])
    err = la.norm(H[0]*C-C*D[0]) + la.norm(H[1]*C-C*D[1]) + la.norm(H[2]*C-C*D[2])
    print(err)
    
    assert CheckMRAPRepresentation(H), "MinimalRepFromMRAP did not return a valid MRAP representation!"
    assert H[0].shape[0]==3 and err<1e-12, "MinimalRepFromMRAP returned a MRAP which is non-similar to the input or has an unexpected size!"
    
    print("H=MinimalRepFromMRAP(D,""obs"")")
    H=MinimalRepFromMRAP(D,"obs")
    print(H[0])
    print(H[1])
    print(H[2])
    
    C = SimilarityMatrix(H[0],D[0])
    err = la.norm(H[0]*C-C*D[0]) + la.norm(H[1]*C-C*D[1]) + la.norm(H[2]*C-C*D[2])
    print(err)
    
    assert CheckMRAPRepresentation(H), "MinimalRepFromMRAP did not return a valid MRAP representation!"
    assert H[0].shape[0]==2 and err<1e-12, "MinimalRepFromMRAP returned a MRAP which is non-similar to the input or has an unexpected size!"
    
    print("H=MinimalRepFromMRAP(D,""obscont"")")
    H=MinimalRepFromMRAP(D,"obscont")
    print(H[0])
    print(H[1])
    print(H[2])
    
    C = SimilarityMatrix(H[0],D[0])
    err = la.norm(H[0]*C-C*D[0]) + la.norm(H[1]*C-C*D[1]) + la.norm(H[2]*C-C*D[2])
    print(err)
    
    assert CheckMRAPRepresentation(H), "MinimalRepFromMRAP did not return a valid MRAP representation!"
    assert H[0].shape[0]==2 and err<1e-12, "MinimalRepFromMRAP returned a MRAP which is non-similar to the input or has an unexpected size!"
    
    print("--SamplesFromMAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-0.17, 0, 0, 0.07],[0.01, -0.78, 0.03, 0.08],[0.22, 0.17, -1.1, 0.02],[0.04, 0.12, 0, -0.42]])
    D1=ml.matrix([[0, 0.06, 0, 0.04],[0.04, 0.19, 0.21, 0.22],[0.22, 0.13, 0.15, 0.19],[0.05, 0, 0.17, 0.04]])
    
    print("Test:")
    print("-----")
    
    print("x=SamplesFromMAP(D0,D1,10000)")
    x=SamplesFromMAP(D0,D1,10000)
    
    print("Moments from the samples:")
    mt = MarginalMomentsFromTrace(x,3)
    print(mt)
    
    print("Moments from the MAP:")
    mm = MarginalMomentsFromMAP(D0,D1,3)
    print(mm)
    
    print("--SamplesFromMMAP--------------------------------------------------------------------------")
    help 
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-1.78, 0.29],[0.07, -0.92]])
    D1=ml.matrix([[0.15, 0.49],[0.23, 0.36]])
    D2=ml.matrix([[0.11, 0.2],[0.01, 0]])
    D3=ml.matrix([[0.14, 0.4],[0.11, 0.14]])
    D = (D0,D1,D2,D3)
    
    print("Test:")
    print("-----")
    
    print("x=SamplesFromMMAP(D,10000)")
    x=SamplesFromMMAP(D,10000)
    
    print("Moments from the samples:")
    mt = MarginalMomentsFromTrace(x[:,0],3)
    print(mt)
    
    print("Moments from the MMAP:")
    mm = MarginalMomentsFromMMAP(D,3)
    print(mm)


if __name__ == "__main__":
    TestMAPPackage()