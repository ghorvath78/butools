# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 17:20:24 2014

@author: gabor
"""

import numpy as np
import numpy.matlib as ml
import butools
from butools.dmap import *
from butools.moments import CheckMoments
from butools.trace import MarginalMomentsFromTrace
from butools.dph import *

def TestDMAPPackage ():

    print("---BuTools: DMAP package test file---")
    
    print("Enable the verbose messages with the BuToolsVerbose flag")
    butools.verbose = True
    
    print("Enable input parameter checking with the BuToolsCheckInput flag")
    butools.checkInput = True
    
    print("--MarginalDistributionFromDMAP--------------------------------------------------------------------")
    help 
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[0, 0.02, 0, 0],[0, 0.17, 0.2, 0.14],[0.16, 0.17, 0.02, 0.18],[0, 0, 0, 0.12]])
    D1=ml.matrix([[0, 0.88, 0.1, 0],[0.18, 0.07, 0.14, 0.1],[0.13, 0.15, 0.15, 0.04],[0.31, 0.18, 0.12, 0.27]])
       
    print("Test:")
    print("-----")
    
    print("[a,A]=MarginalDistributionFromDMAP(D0,D1):")
    [a,A]=MarginalDistributionFromDMAP(D0,D1)
    print(a)
    print(A)
    
    assert a.shape[1]==D0.shape[1] and CheckDPHRepresentation(a,A), "MarginalDistributionFromDMAP returned a wrong DPH representation!"
    
    print("--MarginalMomentsFromDMAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[0, 0.02, 0, 0],[0, 0.17, 0.2, 0.14],[0.16, 0.17, 0.02, 0.18],[0, 0, 0, 0.12]])
    D1=ml.matrix([[0, 0.88, 0.1, 0],[0.18, 0.07, 0.14, 0.1],[0.13, 0.15, 0.15, 0.04],[0.31, 0.18, 0.12, 0.27]])
    
    print("Test:")
    print("-----")
    
    print("moms=MarginalMomentsFromDMAP(D0,D1):")
    moms=MarginalMomentsFromDMAP(D0,D1)
    print(moms)
    
    assert len(moms)==2*D0.shape[1]-1 and CheckMoments(moms), "MarginalMomentsFromDMAP returned wrong moments!"
    
    print("--MarginalDistributionFromDRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[0, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    H1=ml.matrix([[0, 1, -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])

    print("Test:")
    print("-----")
    
    print("[a,A]=MarginalDistributionFromDRAP(H0,H1):")
    [a,A]=MarginalDistributionFromDRAP(H0,H1)
    print(a)
    print(A)
    
    assert a.shape[1]==H0.shape[1] and CheckMGRepresentation(a,A), "MarginalDistributionFromDRAP returned a wrong MG representation!"
    
    print("--MarginalMomentsFromDRAP--------------------------------------------------------------------------")
    help 
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[0, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    H1=ml.matrix([[0, 1, -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])
    
    print("Test:")
    print("-----")
    
    print("moms=MarginalMomentsFromDRAP(H0,H1):")
    moms=MarginalMomentsFromDRAP(H0,H1)
    print(moms)
    
    assert len(moms)==2*H0.shape[1]-1 and CheckMoments(moms), "MarginalMomentsFromDRAP returned wrong moments!"
    
    print("--MarginalDistributionFromDMMAP--------------------------------------------------------------------------")
    help 
    
    D0=ml.matrix([[0.34, 0, 0],[0.06, 0.05, 0.03],[0.11, 0.13, 0]])
    D1=ml.matrix([[0.3, 0, 0],[0.16, 0.18, 0.05],[0.15, 0.04, 0.09]])
    D2=ml.matrix([[0, 0.01, 0],[0.1, 0.07, 0.08],[0.13, 0.12, 0.13]])
    D3=ml.matrix([[0.35, 0, 0],[0, 0.18, 0.04],[0.06, 0.03, 0.01]])

    print("Test:")
    print("-----")
    
    print("[a,A]=MarginalDistributionFromDMMAP((D0,D1,D2,D3)):")
    [a,A]=MarginalDistributionFromDMMAP((D0,D1,D2,D3))
    print(a)
    print(A)
    
    assert a.shape[1]==D0.shape[1] and CheckDPHRepresentation(a,A), "MarginalDistributionFromDMMAP returned a wrong DPH representation!"
    
    print("--MarginalMomentsFromDMMAP--------------------------------------------------------------------------")
    
    D0=ml.matrix([[0.34, 0, 0],[0.06, 0.05, 0.03],[0.11, 0.13, 0]])
    D1=ml.matrix([[0.3, 0, 0],[0.16, 0.18, 0.05],[0.15, 0.04, 0.09]])
    D2=ml.matrix([[0, 0.01, 0],[0.1, 0.07, 0.08],[0.13, 0.12, 0.13]])
    D3=ml.matrix([[0.35, 0, 0],[0, 0.18, 0.04],[0.06, 0.03, 0.01]])
    
    print("Test:")
    print("-----")
    
    print("moms=MarginalMomentsFromDMMAP((D0,D1,D2,D3)):")
    moms=MarginalMomentsFromDMMAP((D0,D1,D2,D3))
    print(moms)
    
    assert len(moms)==2*D0.shape[1]-1 and CheckMoments(moms), "MarginalMomentsFromDMMAP returned wrong moments!"
    
    print("--MarginalDistributionFromDMRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[0.15, 0.2, 0.18],[-0.23, 0.17, 0.22],[0.19, 0.15, 0.16]])
    H1=ml.matrix([[0.01, 0.08, 0.16],[0.02, 0.2, 0.07],[0.02, 0.15, 0.17]])
    H2=ml.matrix([[0.14, 0.07, 0.01],[0.19, 0.02, 0.34],[0.06, 0.1, 0]])

    print("Test:")
    print("-----")
    
    print("[a,A]=MarginalDistributionFromDMRAP((H0,H1,H2)):")
    [a,A]=MarginalDistributionFromDMRAP((H0,H1,H2))
    print(a)
    print(A)
    
    assert a.shape[1]==H0.shape[1] and CheckMGRepresentation(a,A), "MarginalDistributionFromDMRAP returned a wrong MG representation!"
    
    print("--MarginalMomentsFromDMRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[0.15, 0.2, 0.18],[-0.23, 0.17, 0.22],[0.19, 0.15, 0.16]])
    H1=ml.matrix([[0.01, 0.08, 0.16],[0.02, 0.2, 0.07],[0.02, 0.15, 0.17]])
    H2=ml.matrix([[0.14, 0.07, 0.01],[0.19, 0.02, 0.34],[0.06, 0.1, 0]])
    
    print("Test:")
    print("-----")
    
    print("moms=MarginalMomentsFromDMRAP((H0,H1,H2)):")
    moms=MarginalMomentsFromDMRAP((H0,H1,H2))
    print(moms)
    
    assert len(moms)==2*H0.shape[1]-1 and CheckMoments(moms), "MarginalMomentsFromDMRAP returned wrong moments!"
    
    print("--LagCorrelationsFromDMAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[0, 0.02, 0, 0],[0, 0.17, 0.2, 0.14],[0.16, 0.17, 0.02, 0.18],[0, 0, 0, 0.12]])
    D1=ml.matrix([[0, 0.88, 0.1, 0],[0.18, 0.07, 0.14, 0.1],[0.13, 0.15, 0.15, 0.04],[0.31, 0.18, 0.12, 0.27]])
    
    print("Test:")
    print("-----")
    
    print("LagCorrelationsFromDMAP(D0,D1,3):")
    corr = LagCorrelationsFromDMAP(D0,D1,3)
    print(corr)
    
    assert len(corr)==3 and np.all(np.array(corr)<1) and np.all(np.array(corr)>-1), "LagCorrelationsFromDMAP returned wrong autocorrelation coefficients!"
    
    print("--LagCorrelationsFromDRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[0, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    H1=ml.matrix([[0, 1, -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])
    
    print("Test:")
    print("-----")
    
    print("LagCorrelationsFromDRAP(H0,H1,3):")
    corr = LagCorrelationsFromDRAP(H0,H1,3)
    print(corr)
    
    assert len(corr)==3 and np.all(np.array(corr)<1) and np.all(np.array(corr)>-1), "LagCorrelationsFromDRAP returned wrong autocorrelation coefficients!"
    
    print("--LagkJointMomentsFromDMAP--------------------------------------------------------------------------")
    help 
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[0, 0.02, 0, 0],[0, 0.17, 0.2, 0.14],[0.16, 0.17, 0.02, 0.18],[0, 0, 0, 0.12]])
    D1=ml.matrix([[0, 0.88, 0.1, 0],[0.18, 0.07, 0.14, 0.1],[0.13, 0.15, 0.15, 0.04],[0.31, 0.18, 0.12, 0.27]])
    
    print("Test:")
    print("-----")
    
    print("LagkJointMomentsFromDMAP(D0,D1,2,1):")
    Nm=LagkJointMomentsFromDMAP(D0,D1,4,1)
    print(Nm)
    
    moms=MarginalMomentsFromDMAP(D0,D1,4)
    print(moms)
    
    print("Correlation from joint moments:")
    cjm=np.empty(3)
    for i in range(3):
        Nx=LagkJointMomentsFromDMAP(D0,D1,1,i+1)
        cjm[i] = (Nx[1,1]-moms[0]**2) / (moms[1]-moms[0]**2)
    print(cjm)
    corr = LagCorrelationsFromDMAP(D0,D1,3)
    print(corr)    
    
    assert np.all(np.all(Nm>0)) and la.norm(moms-Nm[0,1:])<1e-14 and la.norm(moms-Nm[1:,0].A.flatten())<1e-14 and la.norm(corr-cjm)<1e-14, "Joint moment matrix is invalid!"
    

    print("--LagkJointMomentsFromDRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[0, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    H1=ml.matrix([[0, 1, -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])
    
    print("Test:")
    print("-----")
    
    print("LagkJointMomentsFromDRAP(H0,H1,2,1):")
    Nm=LagkJointMomentsFromDRAP(H0,H1,4,1)
    print(len(Nm))
    
    
    moms=MarginalMomentsFromDRAP(H0,H1,4)
    
    print("Correlation from joint moments:")
    cjm=np.empty(3)
    for i in range(3):
        Nx=LagkJointMomentsFromDRAP(H0,H1,1,i+1)
        cjm[i] = (Nx[1,1]-moms[0]**2) / (moms[1]-moms[0]**2)
    print(cjm)
    corr = LagCorrelationsFromDRAP(H0,H1,3)
    
    assert np.all(np.all(Nm>0)) and la.norm(moms-Nm[0,1:])<1e-14 and la.norm(moms-Nm[1:,0].A.flatten())<1e-14 and la.norm(corr-cjm)<1e-14, "Joint moment matrix is invalid!"
    
    print("--LagkJointMomentsFromDMMAP--------------------------------------------------------------------------")
    help 
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[0.34, 0, 0],[0.06, 0.05, 0.03],[0.11, 0.13, 0]])
    D1=ml.matrix([[0.3, 0, 0],[0.16, 0.18, 0.05],[0.15, 0.04, 0.09]])
    D2=ml.matrix([[0, 0.01, 0],[0.1, 0.07, 0.08],[0.13, 0.12, 0.13]])
    D3=ml.matrix([[0.35, 0, 0],[0, 0.18, 0.04],[0.06, 0.03, 0.01]])
    
    print("Test:")
    print("-----")
    
    print("LagkJointMomentsFromDMMAP((D0,D1,D2,D3),3,1)")
    Nm=LagkJointMomentsFromDMMAP((D0,D1,D2,D3),3,1)
    
    print("Moments of arrival type 1, Nm[0]:")
    print(Nm[0])
    print("Moments of arrival type 2, Nm[1]:")
    print(Nm[1])
    print("Moments of arrival type 2, Nm[2]:")
    print(Nm[2])
    
    assert len(Nm)==3 and la.norm(SumMatrixList(Nm)-LagkJointMomentsFromDMAP(D0,D1+D2+D3,3,1))<1e-13,"Joint moment matrix is invalid!"
    
    print("--LagkJointMomentsFromDMRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[0.15, 0.2, 0.18],[-0.23, 0.17, 0.22],[0.19, 0.15, 0.16]])
    H1=ml.matrix([[0.01, 0.08, 0.16],[0.02, 0.2, 0.07],[0.02, 0.15, 0.17]])
    H2=ml.matrix([[0.14, 0.07, 0.01],[0.19, 0.02, 0.34],[0.06, 0.1, 0]])

    print("Test:")
    print("-----")
    
    print("LagkJointMomentsFromDMRAP((H0,H1,H2),3,2)")
    Nm=LagkJointMomentsFromDMRAP((H0,H1,H2),3,2)
    
    print("Moments of arrival type 1, Nm[0]:")
    print(Nm[0])
    print("Moments of arrival type 2, Nm[1]:")
    print(Nm[1])
    
    assert len(Nm)==2 and la.norm(SumMatrixList(Nm)-LagkJointMomentsFromDRAP(H0,H1+H2,3,2))<1e-13,"Joint moment matrix is invalid!"
    
    print("--RandomDMAP--------------------------------------------------------------------------")
    
    print("Test:")
    print("-----")
    
    print("[D0,D1]=RandomDMAP(4,5.62,10):")
    [D0,D1]=RandomDMAP(4,5.62,10)
    print(D0)
    print(D1)
    
    print("Check the mean of the obtained DMAP:")
    print(np.sum(D0+D1,1)-1.0)
    m = MarginalMomentsFromDMAP(D0,D1,1)[0]
    print(m)
    
    assert CheckDMAPRepresentation(D0,D1), "RandomDMAP failed to return a valid DMAP representation!"
    assert np.max(np.abs(m-5.62))<1e-14, "RandomDMAP failed to match the given mean value!"
    
    print("--RandomDMMAP--------------------------------------------------------------------------")
    
    print("Test:")
    print("-----")
    
    print("D=RandomDMMAP(4,3,5.62,10):")
    D=RandomDMMAP(4,3,5.62,10)
    print(D[0])
    print(D[1])
    print(D[2])
    print(D[3])
    
    print("Check the mean of the obtained DMMAP:")
    m = MarginalMomentsFromDMMAP(D,1)[0]
    print(m)
    
    assert CheckDMMAPRepresentation(D), "RandomDMMAP failed to return a valid DMMAP representation!"
    assert np.max(np.abs(m-5.62))<1e-14, "RandomDMMAP failed to match the given mean value!"
    
    print("-CheckDMAPRepresentation---------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[0, 0.02, 0],[0, 0.17, 0.2],[0.16, 0.17, 0.02]])
    D1=ml.matrix([[0, 0.88, 0.1, 0],[0.18, 0.07, 0.14, 0.1],[0.13, 0.15, 0.15, 0.04],[0.31, 0.18, 0.12, 0.27]])
    
    print("Test:")
    print("-----")
    
    print("CheckDMAPRepresentation(D0,D1):")
    flag=CheckDMAPRepresentation(D0,D1)
    print(flag)
    
    assert flag==False, "CheckDMAPRepresentation failed to detect the incompatible shapes of D0 and D1!"
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[0, 0.02, 0],[0, 0.17, 0.2],[0.16, 0.17, 0.02]])
    D1=ml.matrix([[0, 0.88, 0.1],[0.18, 0.07, 0.14],[0.13, 0.15, 0.15]])
    
    print("Test:")
    print("-----")
    
    print("CheckDMAPRepresentation(D0,D1):")
    flag=CheckDMAPRepresentation(D0,D1)
    print(flag)
    
    assert flag==False, "CheckDMAPRepresentation failed to detect invalid rowsums!"
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[0, 0.02, 0, 0],[0, 0.17, 0.2, 0.14],[0.16, 0.17, 0.02, 0.18],[0, 0, 0, 0.12]])
    D1=ml.matrix([[0, 0.88, 0.1, 0],[0.18, 0.07, 0.14, 0.1],[0.13, 0.15, 0.15, 0.04],[0.31, 0.18, 0.12, 0.27]])
    
    print("Test:")
    print("-----")
    
    print("CheckDMAPRepresentation(D0,D1):")
    flag=CheckDMAPRepresentation(D0,D1)
    print(flag)
    
    assert flag==True, "CheckDMAPRepresentation failed to recognize a valid DMAP representation!"
    
    print("--CheckDRAPRepresentation--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[0, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02],[0.2, 0, 0]])
    H1=ml.matrix([[0, 1, -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29],[0, 0.8, 0]])
    
    print("Test:")
    print("-----")
    
    print("CheckDRAPRepresentation(H0,H1):")
    flag=CheckDRAPRepresentation(H0,H1)
    print(flag)
    
    assert flag==False, "CheckDRAPRepresentation failed to detect the incompatible shapes of D0 and D1!"
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[0.2, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    H1=ml.matrix([[0, 1, -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])
    
    print("Test:")
    print("-----")
    
    print("CheckDRAPRepresentation(H0,H1):")
    flag=CheckDRAPRepresentation(H0,H1)
    print(flag)
    
    assert flag==False, "CheckDRAPRepresentation failed to detect invalid rowsums!"
    
    print("Input:")
    print("------")
    
    x=15    
    H0=ml.matrix([[0, 0, x],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    H1=ml.matrix([[0, 1, -x],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])
    
    print("Test:")
    print("-----")
    
    print("CheckDRAPRepresentation(H0,H1):")
    flag=CheckDRAPRepresentation(H0,H1)
    print(flag)
    
    assert flag==False, "CheckDRAPRepresentation failed to detect invalid eigenvalues!"
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[0, 0.5, 0.1],[0, -1.4, 3.1],[0.67, 0, 0.4]])
    H1=ml.matrix([[0, 0.4, 0],[0, -0,2, -0.5],[0.3, -0.7, 0.33]])

    print("Test:")
    print("-----")
    
    print("CheckDRAPRepresentation(H0,H1):")
    flag=CheckDRAPRepresentation(H0,H1)
    print(flag)
    
    assert flag==False, "CheckDRAPRepresentation failed to detect invalid eigenvalues!"
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[0, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    H1=ml.matrix([[0, 1, -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])
    
    print("Test:")
    print("-----")
    
    print("CheckDRAPRepresentation(H0,H1):")
    flag=CheckDRAPRepresentation(H0,H1)
    print(flag)
    
    assert flag==True, "CheckDRAPRepresentation failed to recognize a valid DRAP representation!"
    
    print("-CheckDMMAPRepresentation---------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[0.34, 0, 0],[0.06, 0.05, 0.03],[0.11, 0.13, 0]])
    D1=ml.matrix([[0.3, 0, 0],[0.16, 0.18, 0.05],[0.15, 0.04, 0.09]])
    D2=ml.matrix([[0, 0.01, 0],[0.1, 0.07, 0.08],[0.13, 0.12, 0.13]])
    D3=ml.matrix([[0.35, 0, 0],[0, 0.18, 0.04],[0.06, 0.03, 0.01]])
    
    print("Test:")
    print("-----")
    
    print("CheckDMMAPRepresentation((D0,D1,D2,D3)):")
    flag=CheckDMMAPRepresentation((D0,D1,D2,D3))
    print(flag)
    
    assert flag==True, "CheckDMMAPRepresentation failed to recognize a valid DMMAP representation!"
    
    print("-CheckDMRAPRepresentation---------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[0.15, 0.2, 0.18],[-0.23, 0.17, 0.22],[0.19, 0.15, 0.16]])
    H1=ml.matrix([[0.01, 0.08, 0.16],[0.02, 0.2, 0.07],[0.02, 0.15, 0.17]])
    H2=ml.matrix([[0.14, 0.07, 0.01],[0.19, 0.02, 0.34],[0.06, 0.1, 0]])
    
    print("Test:")
    print("-----")
    
    print("CheckDMRAPRepresentation((H0,H1,H2)):")
    flag=CheckDMRAPRepresentation((H0,H1,H2))
    print(flag)
    
    assert flag==True, "CheckDMRAPRepresentation failed to recognize a valid DMRAP representation!"
    
    print("--DRAPFromMoments--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    G0=ml.matrix([[0, 0.02, 0, 0],[0, 0.17, 0.2, 0.14],[0.16, 0.17, 0.02, 0.18],[0, 0, 0, 0.12]])
    G1=ml.matrix([[0, 0.88, 0.1, 0],[0.18, 0.07, 0.14, 0.1],[0.13, 0.15, 0.15, 0.04],[0.31, 0.18, 0.12, 0.27]])
    moms=MarginalMomentsFromDRAP(G0,G1,5)
    Nm=LagkJointMomentsFromDRAP(G0,G1,2,1)
    print(moms)
    print(Nm)
    
    print("Test:")
    print("-----")
    
    print("[H0,H1]=DRAPFromMoments(moms,Nm):")
    [H0,H1]=DRAPFromMoments(moms,Nm)
    print(H0)
    print(H1)
    
    rmoms=MarginalMomentsFromDRAP(H0,H1,5,1e-12)
    rNm=LagkJointMomentsFromDRAP(H0,H1,2,1,1e-12)
    print(rmoms)
    print(rNm)
    
    assert la.norm(np.array(moms)-np.array(rmoms))<1e-11 and la.norm(Nm-rNm)<1e-12, "The moments and joint moments returned by DRAPFromMoments are not the same as given!"
    
    print("--DMRAPFromMoments--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")

    G0=ml.matrix([[0.34, 0, 0],[0.06, 0.05, 0.03],[0.11, 0.13, 0]])
    G1=ml.matrix([[0.3, 0, 0],[0.16, 0.18, 0.05],[0.15, 0.04, 0.09]])
    G2=ml.matrix([[0, 0.01, 0],[0.1, 0.07, 0.08],[0.13, 0.12, 0.13]])
    G3=ml.matrix([[0.35, 0, 0],[0, 0.18, 0.04],[0.06, 0.03, 0.01]])
    G=(G0,G1,G2,G3)
    moms=MarginalMomentsFromDMRAP(G,5)
    print(moms)
    Nm=LagkJointMomentsFromDMRAP(G,2,1)
    Nm1, Nm2, Nm3 = Nm
    print(Nm1)
    print(Nm2)
    print(Nm3)
    
    print("Test:")
    print("-----")
    
    print("H=DMRAPFromMoments(moms,Nm):")
    H=DMRAPFromMoments(moms,Nm)
    print(H[0])
    print(H[1])
    print(H[2])
    print(H[3])
    
    rmoms=MarginalMomentsFromDMRAP(H,5,1e-10)
    print(rmoms)
    rNm1, rNm2, rNm3 = LagkJointMomentsFromDMRAP(H,2,1,1e-10)
    print(rNm1)    
    print(rNm2)    
    print(rNm3)    
    
    assert la.norm(np.array(moms)-np.array(rmoms))<1e-8 and la.norm(Nm1-rNm1)<1e-9 and la.norm(Nm2-rNm2)<1e-9 and la.norm(Nm3-rNm3)<1e-9, "The moments and joint moments returned by DMRAPFromMoments are not the same as given!"
    
    print("--DMAPFromDRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[0, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    H1=ml.matrix([[0, 1, -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])
   
    print("Test:")
    print("-----")
    
    print("[H0,H1]=DMAPFromDRAP(H0,H1):")
    [D0,D1]=DMAPFromDRAP(H0,H1)
    print(D0)
    print(D1)
    
    err = la.norm(LagkJointMomentsFromDRAP(D0,D1,3,1)-LagkJointMomentsFromDRAP(H0,H1,3,1))
    print(err)
    assert err<1e-10, "The DMAP returned by DMAPFromDRAP is not similar to the input!"
    assert CheckDMAPRepresentation(D0,D1), "The result of DMAPFromDRAP is not a DMAP, as it should be!"
    
    print("--DMMAPFromDMRAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    H0=ml.matrix([[0.15, 0.2, 0.18],[-0.23, 0.17, 0.22],[0.19, 0.15, 0.16]])
    H1=ml.matrix([[0.01, 0.08, 0.16],[0.02, 0.2, 0.07],[0.02, 0.15, 0.17]])
    H2=ml.matrix([[0.14, 0.07, 0.01],[0.19, 0.02, 0.34],[0.06, 0.1, 0]])
    H=(H0,H1,H2)
    
    moms=MarginalMomentsFromDMRAP(H)
    jmom=LagkJointMomentsFromDMRAP(H,3,1)
    print(moms)   
    print(jmom)
    
    print("Test:")
    print("-----")
    
    print("G=DMMAPFromDMRAP({H0,H1,H2}):")
    G=DMMAPFromDMRAP(H)
    print(G[0])
    print(G[1])
    print(G[2])
    
    rmoms=MarginalMomentsFromDMMAP(G)
    rjmom=LagkJointMomentsFromDMMAP(G,3,1)
    print(rmoms)   
    print(rjmom)
    
    err = la.norm(rjmom[0]-jmom[0]) + la.norm(rjmom[1]-jmom[1])
    print(err)
    assert err<1e-12, "The DMMAP returned by DMMAPFromDMRAP is not similar to the input!"
    assert CheckDMMAPRepresentation(G), "The result of DMMAPFromDMRAP is not a DMMAP, as it should be!"
    
    print("--CanonicalFromDMAP2--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[0.46, 0.28],[0.35, 0.23]])
    D1=ml.matrix([[0.08, 0.18],[0.14, 0.28]])
    
    print("Test:")
    print("-----")
    
    print("[H0,H1]=CanonicalFromDMAP2(D0,D1):")
    [H0,H1]=CanonicalFromDMAP2(D0,D1)
    print(H0)
    print(H1)
    
    C=SimilarityMatrix(H0,D0)
    err = la.norm(H0*C-C*D0) + la.norm(H1*C-C*D1)
    print(err)
    
    assert CheckDMAPRepresentation(H0,H1), "The result of CanonicalFromDMAP2 is not a valid DMAP representation!"
    assert err<1e-12, "The MAP returned by CanonicalFromDMAP2 is not similar to the input!"      
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[0.26, 0.28],[0.35, 0.23]])
    D1=ml.matrix([[0.28, 0.18],[0.14, 0.28]])
    
    print("Test:")
    print("-----")
    
    print("[H0,H1]=CanonicalFromDMAP2(D0,D1):")
    [H0,H1]=CanonicalFromDMAP2(D0,D1)
    print(H0)
    print(H1)
    
    C=SimilarityMatrix(H0,D0)
    err = la.norm(H0*C-C*D0) + la.norm(H1*C-C*D1)
    print(err)
    
    assert CheckDMAPRepresentation(H0,H1), "The result of CanonicalFromDMAP2 is not a valid DMAP representation!"
    assert err<1e-12, "The MAP returned by CanonicalFromDMAP2 is not similar to the input!"      
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[0.14, 0.34],[0.35, 0.23]])
    D1=ml.matrix([[0.22, 0.3],[0.28, 0.14]])
    
    print("Test:")
    print("-----")
    
    print("[H0,H1]=CanonicalFromDMAP2(D0,D1):")
    [H0,H1]=CanonicalFromDMAP2(D0,D1)
    print(H0)
    print(H1)
    
    C=SimilarityMatrix(H0,D0)
    err = la.norm(H0*C-C*D0) + la.norm(H1*C-C*D1)
    print(err)
    
    assert CheckDMAPRepresentation(H0,H1), "The result of CanonicalFromDMAP2 is not a valid DMAP representation!"
    assert err<1e-12, "The MAP returned by CanonicalFromDMAP2 is not similar to the input!"      
    
    print("--DMAP2FromMoments--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[0.2, 0.7],[0.6, 0.1]])
    D1=ml.matrix([[0.09, 0.01],[0.2, 0.1]])
    moms=MarginalMomentsFromDMAP(D0, D1, 3)
    corr=LagCorrelationsFromDMAP(D0, D1, 1)
    print(moms)
    print(corr)
    
    print("Test:")
    print("-----")
    
    print("[D0,D1]=DMAP2FromMoments(moms,corr):")
    [D0,D1]=DMAP2FromMoments(moms,corr)
    print(D0)
    print(D1)
    
    rmoms=MarginalMomentsFromDMAP(D0, D1, 3)
    rcorr=LagCorrelationsFromDMAP(D0, D1, 1)
    print(rmoms)
    print(rcorr)    
    
    assert CheckDMAPRepresentation(D0,D1), "DMAP2FromMoments returned an invalid DMAP representation!"
    assert la.norm(np.array(moms)-np.array(rmoms))<1e-11 and abs(corr-rcorr)<1e-11, "The moments and the correlation returned by DMAP2FromMoments are not the same as given!"

    print("--SamplesFromDMAP--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[0, 0.02, 0, 0],[0, 0.17, 0.2, 0.14],[0.16, 0.17, 0.02, 0.18],[0, 0, 0, 0.12]])
    D1=ml.matrix([[0, 0.88, 0.1, 0],[0.18, 0.07, 0.14, 0.1],[0.13, 0.15, 0.15, 0.04],[0.31, 0.18, 0.12, 0.27]])

    print("Test:")
    print("-----")
    
    print("x=SamplesFromDMAP(D0,D1,10000)")
    x=SamplesFromDMAP(D0,D1,10000)
    
    print("Moments from the samples:")
    mt = MarginalMomentsFromTrace(x,3)
    print(mt)
    
    print("Moments from the DMAP:")
    mm = MarginalMomentsFromDMAP(D0,D1,3)
    print(mm)
    
    print("--SamplesFromDMMAP--------------------------------------------------------------------------")
    help 
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[0.34, 0, 0],[0.06, 0.05, 0.03],[0.11, 0.13, 0]])
    D1=ml.matrix([[0.3, 0, 0],[0.16, 0.18, 0.05],[0.15, 0.04, 0.09]])
    D2=ml.matrix([[0, 0.01, 0],[0.1, 0.07, 0.08],[0.13, 0.12, 0.13]])
    D3=ml.matrix([[0.35, 0, 0],[0, 0.18, 0.04],[0.06, 0.03, 0.01]])
    D = (D0,D1,D2,D3)
    
    print("Test:")
    print("-----")
    
    print("x=SamplesFromDMMAP(D,10000)")
    x=SamplesFromDMMAP(D,10000)
    
    print("Moments from the samples:")
    mt = MarginalMomentsFromTrace(x[:,0],3)
    print(mt)
    
    print("Moments from the DMMAP:")
    mm = MarginalMomentsFromDMMAP(D,3)
    print(mm)


if __name__ == "__main__":
    TestDMAPPackage()