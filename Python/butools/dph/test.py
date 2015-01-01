# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 17:20:24 2014

@author: gabor
"""

import numpy as np
import numpy.matlib as ml
import butools
from butools.dph import *
from butools.moments import CheckMoments
from butools.trace import MarginalMomentsFromTrace

def TestDPHPackage ():

    print('---BuTools: DPH package test file---')
    
    print('Enable the verbose messages with the BuToolsVerbose flag')
    butools.verbose = True
    
    print('Enable input parameter checking with the BuToolsCheckInput flag')
    butools.checkInput = True
    
    print("--MomentsFromMG---------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a = ml.matrix([[-0.6, 0.3, 1.3]])
    print(a)
    A = ml.matrix([[0.25, 0.2, -0.15],[0.3, 0.1, 0.25],[0, 0.2, 0.47]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("MomentsFromMG(a,A,3):")
    moms=MomentsFromMG(a,A)
    print(moms)
    
    assert np.all(np.array(moms)>0) and CheckMoments(moms)==True, "The function returned invalid moments!"
    
    print("Test:")
    print("-----")
    
    print("MomentsFromMG(a,A):")
    moms=MomentsFromMG(a,A,3)
    print(moms)
    
    assert np.all(np.array(moms)>0) and CheckMoments(moms)==True, "The function returned invalid moments!"
    
    print("--MomentsFromDPH---------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.76, 0, 0.24]])
    print(a)
    A=ml.matrix([[0.34, 0.66, 0],[0.79, 0.05, 0.07],[0.26, 0.73, 0.01]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("MomentsFromDPH(a,A,5)")
    moms=MomentsFromDPH(a,A,5)
    print(moms)
    
    assert np.all(np.array(moms)>0) and CheckMoments(moms)==True, "The function returned invalid moments!"

    print("--PmfFromMG-------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a = ml.matrix([[-0.6, 0.3, 1.3]])
    print(a)
    A = ml.matrix([[0.25, 0.2, -0.15],[0.3, 0.1, 0.25],[0, 0.2, 0.47]])
    print(A)
    x = np.linspace(0,100,101);
    
    print("Test:")
    print("-----")
    
    print("pmf=PmfFromMG(a,A,x):")
    pmf = PmfFromMG(a,A,x)
    
    assert np.all(pmf)>=0, "PmfFromMG returned negative pmf!"
    assert abs(np.dot(pmf,x) - MomentsFromMG(a,A,1)[0])<1e-12, "The mean computed from the pmf does not match the theoretical result!"
    
    print("--PmfFromDPH-------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.76, 0, 0.24]])
    print(a)
    A=ml.matrix([[0.34, 0.66, 0],[0.79, 0.05, 0.07],[0.26, 0.73, 0.01]])
    print(A)
    x = np.linspace(0,1000,1001)
    
    print("Test:")
    print("-----")
    
    print("pmf=PmfFromDPH(a,A,x):")
    pmf = PmfFromDPH(a,A,x)
    
    assert np.all(pmf)>=0, "PmfFromDPH returned negative pmf!"
    assert abs(np.dot(pmf,x) - MomentsFromDPH(a,A,1))<1e-12, "The mean computed from the pmf does not match the theoretical result!"
    
    print("--CdfFromMG-------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a = ml.matrix([[-0.6, 0.3, 1.3]])
    print(a)
    A = ml.matrix([[0.25, 0.2, -0.15],[0.3, 0.1, 0.25],[0, 0.2, 0.47]])
    print(A)
    x = np.linspace(0,100,101);
    
    print("Test:")
    print("-----")
    
    print("cdf=CdfFromMG(a,A,x):")
    cdf = CdfFromMG(a,A,x)
    
    assert np.all(np.diff(cdf)>=0), "The cdf is not increasing monotonously!"
    assert abs(np.sum(1.0-cdf) - MomentsFromMG(a,A,1))<1e-12, "The mean computed from the cdf does not match the theoretical result!"
    
    print("--CdfFromDPH-------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.76, 0, 0.24]])
    print(a)
    A=ml.matrix([[0.34, 0.66, 0],[0.79, 0.05, 0.07],[0.26, 0.73, 0.01]])
    print(A)
    x = np.linspace(0,1000,1001)
    
    print("Test:")
    print("-----")
    
    print("cdf=CdfFromDPH(a,A,x):")
    cdf = CdfFromDPH(a,A,x)
    
    assert np.all(np.diff(cdf)>=0), "The cdf is not increasing monotonously!"
    assert abs(np.sum(1.0-cdf) - MomentsFromDPH(a,A,1))<1e-12, "The mean computed from the cdf does not match the theoretical result!"
       
    print("--RandomDPH--------------------------------------------------------------------------")
    
    print("Test:")
    print("-----")
       
    print("[a,A]=RandomDPH(3,10,5):")
    a,A=RandomDPH(3,10,5)
    print(a)
    print(A)
       
    assert CheckDPHRepresentation(a,A), "RandomDPH failed to return a valid DPH representation!"
    assert np.max(np.abs(MomentsFromDPH(a,A,1)[0]-10.0))<1e-14, "RandomDPH failed to match the given mean value!"
    clo = 1-np.sum(A,1)
    assert np.sum(a==0)+np.sum(A==0)+np.sum(clo==0)==5, "The number of zero entries does not match the function parameter!"
    
    print("--CheckMGRepresentation-------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a = ml.matrix([[-0.6, 0.3, 1.3]])
    print(a)
    A = ml.matrix([[0.25, 0.2, -0.15],[0.3, 0.1, 0.25],[0, 0.2, 0.47]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("CheckMGRepresentation(a,A):")
    flag=CheckMGRepresentation(a,A)
    print(flag)
    
    assert flag==True, "CheckMGRepresentation failed to recognize a valid MG distribution!"
    
    print("Input:")
    print("------")
    
    a = ml.matrix([[-0.6, 0.3, 1.3]])
    print(a)
    A = ml.matrix([[0.35, 0.2, -0.25],[0.3, 0.1, 0.25],[0, 0.2, 0.47]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("CheckMGRepresentation(a,A):")
    flag=CheckMGRepresentation(a,A)
    print(flag)
    
    assert flag==False, "CheckMGRepresentation failed to recognize wrong eigenvalues!"
       
    print("--CheckDPHRepresentation-------------------------------------------------------------")
    
    print("Input:")
    print("------")

    a=ml.matrix([[0.48, 0.08, 0.26, 0.18]])
    print(a)
    A=ml.matrix([[0, 0.08, 0.08, 0.8],[0.55, 0, 0.24, 0.19],[0.06, 0.03, 0, 0.001],[0.23, 0.005, 0.2, 0.53]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("CheckDPHRepresentation(a,A):")
    flag=CheckDPHRepresentation(a,A)
    print(flag)
    
    assert flag==True, "CheckDPHRepresentation failed to recognize a valid DPH distribution!"
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.48, 0.08]])
    print(a)
    A=ml.matrix([[0, 0.08],[0.55, 0.5]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("CheckDPHRepresentation(a,A):")
    flag=CheckDPHRepresentation(a,A)
    print(flag)
    
    assert flag==False, "CheckDPHRepresentation failed to recognize wrong row sums!"
    
    print("--MGFromMoments---------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    moms = [4.08, 20.41, 130.45, 1054.41, 10463.73]
    print(moms)
    
    print("Test:")
    print("-----")
    
    print("[a,A]=MGFromMoments(moms):")
    a,A=MGFromMoments(moms)
    
    print("MomentsFromMG(a,A,5):")
    memoms = MomentsFromMG(a,A,5)
    print(memoms)
    
    assert la.norm(np.array(moms)-np.array(memoms))<1e-9, "The moments of the result of MGFromMoments do not match the input!"
    
    print("--DPHFromMG--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[-0.6, 0.3, 1.3]])
    print(a)
    A=ml.matrix([[0.1, 0.2, 0],[0.3, 0.1, 0.25],[-0.3, 0.2, 0.77]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("CheckMGRepresentation(a,A):")
    flag=CheckMGRepresentation(a,A)
    print(flag)
    
    print("CheckDPHRepresentation(a,A):")
    flag=CheckDPHRepresentation(a,A)
    print(flag)
    
    print("[b,B]=DPHFromMG(a,A)")
    b,B=DPHFromMG(a,A)
    print(b)
    print(B)
    
    print("Check the obtained DPH, CheckDPHRepresentation(b,B):")
    flag=CheckDPHRepresentation(b,B)
    print(flag)
    
    C=SimilarityMatrix(A,B)
    err1 = np.max(np.max(np.abs(A*C-C*B)))
    err2 = np.max(np.abs(a*C-b))
    print("Transformation errors:")
    print(err1)
    print(err2)
    
    assert flag and err1<1e-12 and err2<1e-12, "Transformation to DPH failed!"
    
    print("--CanonicalFromDPH2------------------------------------------------------------------")
    
    print("Input:")
    print("------")
       
    a=ml.matrix([[0, 1.0]])
    A=ml.matrix([[0.23, 0.22],[0.41, 0.48]])
    
    print("Test:")
    print("-----")
    
    print("[b,B]=CanonicalFromDPH2(a,A):")
    b,B=CanonicalFromDPH2(a,A)
    print(b)
    print(B)
    
    print("Eigenvalues of A:")
    ev = la.eigvals(A)
    ix = np.argsort(np.abs(np.real(ev)))
    print(ev[ix])
    
    print("Check the obtained DPH, CheckDPHRepresentation(b,B):")
    flag=CheckDPHRepresentation(b,B)
    print(flag)   
    
    C=SimilarityMatrix(A,B)
    err1 = np.max(np.max(np.abs(A*C-C*B)))
    err2 = np.max(np.abs(a*C-b))
    
    print("Transformation errors:")
    print(err1)
    print(err2)
    
    assert err1<1e-12 and err2<1e-12, "Transformation to canonical DPH(2) failed!"
    
    print("Input:")
    print("------")
       
    a=ml.matrix([[1.0, 0]])
    A=ml.matrix([[0, 0.61],[0.56, 0.44]])
    
    print("Test:")
    print("-----")
    
    print("[b,B]=CanonicalFromDPH2(a,A):")
    b,B=CanonicalFromDPH2(a,A)
    print(b)
    print(B)
    
    print("Eigenvalues of A:")
    ev = la.eigvals(A)
    ix = np.argsort(np.abs(np.real(ev)))
    print(ev[ix])
    
    print("Check the obtained DPH, CheckDPHRepresentation(b,B):")
    flag=CheckDPHRepresentation(b,B)
    print(flag)   
    
    C=SimilarityMatrix(A,B)
    err1 = np.max(np.max(np.abs(A*C-C*B)))
    err2 = np.max(np.abs(a*C-b))
    
    print("Transformation errors:")
    print(err1)
    print(err2)
    
    assert err1<1e-12 and err2<1e-12, "Transformation to canonical DPH(2) failed!"

    print("--CanonicalFromDPH3------------------------------------------------------------------")
    
    print("Input:")
    print("------")

    a=ml.matrix([[0.46, 0.22, 0.32]])
    print(a)
    A=ml.matrix([[0.67, 0.01, 0.12],[0.06, 0.45, 0.15],[0.18, 0.43, 0.32]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("[b,B]=CanonicalFromDPH3(a,A):")
    b,B=CanonicalFromDPH3(a,A)
    print(b)
    print(B)
   
    print("Eigenvalues of A:")
    ev = la.eigvals(A)
    ix = np.argsort(np.abs(np.real(ev)))
    print(ev[ix])
    
    print("Check the obtained DPH, CheckDPHRepresentation(b,B):")
    flag=CheckDPHRepresentation(b,B)
    print(flag)     
   
    C=SimilarityMatrix(A,B)
    err1 = np.max(np.max(np.abs(A*C-C*B)))
    err2 = np.max(np.abs(a*C-b))
    
    print("Transformation errors:")
    print(err1)
    print(err2)
    
    assert err1<1e-12 and err2<1e-12, "Transformation to canonical PH(3) failed!"
 
    print("Input:")
    print("------")

    a=ml.matrix([[0.76, 0.12, 0.12]])
    print(a)
    A=ml.matrix([[0.31, 0., 0.],[0.98, 0., 0.02],[0.88, 0.04, 0.08]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("[b,B]=CanonicalFromDPH3(a,A):")
    b,B=CanonicalFromDPH3(a,A)
    print(b)
    print(B)

    print("Eigenvalues of A:")
    ev = la.eigvals(A)
    ix = np.argsort(np.abs(np.real(ev)))
    print(ev[ix])
    
    print("Check the obtained DPH, CheckDPHRepresentation(b,B):")
    flag=CheckDPHRepresentation(b,B)
    print(flag)     
   
    C=SimilarityMatrix(A,B)
    err1 = np.max(np.max(np.abs(A*C-C*B)))
    err2 = np.max(np.abs(a*C-b))
    
    print("Transformation errors:")
    print(err1)
    print(err2)
    
    assert err1<1e-12 and err2<1e-12, "Transformation to canonical PH(3) failed!"
 
    print("Input:")
    print("------")

    a=ml.matrix([[0.67, 0.07, 0.26]])
    print(a)
    A=ml.matrix([[0.31, 0., 0.],[0.98, 0., 0.02],[0.88, 0.04, 0.08]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("[b,B]=CanonicalFromDPH3(a,A):")
    b,B=CanonicalFromDPH3(a,A)
    print(b)
    print(B)
   
    print("Eigenvalues of A:")
    ev = la.eigvals(A)
    ix = np.argsort(np.abs(np.real(ev)))
    print(ev[ix])
    
    print("Check the obtained DPH, CheckDPHRepresentation(b,B):")
    flag=CheckDPHRepresentation(b,B)
    print(flag)     
   
    C=SimilarityMatrix(A,B)
    err1 = np.max(np.max(np.abs(A*C-C*B)))
    err2 = np.max(np.abs(a*C-b))
    
    print("Transformation errors:")
    print(err1)
    print(err2)
    
    assert err1<1e-12 and err2<1e-12, "Transformation to canonical PH(3) failed!"
 
    print("Input:")
    print("------")

    a=ml.matrix([[0.78, 0.04, 0.18]])
    print(a)
    A=ml.matrix([[0.06, 0.25, 0.31],[0.45, 0.18, 0.33],[0.98, 0, 0.01]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("[b,B]=CanonicalFromDPH3(a,A):")
    b,B=CanonicalFromDPH3(a,A)
   
    print("Eigenvalues of A:")
    ev = la.eigvals(A)
    ix = np.argsort(np.abs(np.real(ev)))
    print(ev[ix])
    
    print("Check the obtained DPH, CheckDPHRepresentation(b,B):")
    flag=CheckDPHRepresentation(b,B)
    print(flag)     
   
    C=SimilarityMatrix(A,B)
    err1 = np.max(np.max(np.abs(A*C-C*B)))
    err2 = np.max(np.abs(a*C-b))
    
    print("Transformation errors:")
    print(err1)
    print(err2)
    
    assert err1<1e-12 and err2<1e-12, "Transformation to canonical PH(3) failed!"

    print("--AcyclicDPHFromMG-------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0, 0, 1.0]])
    A=ml.matrix([[0.22, 0, 0],[0.3, 0.1, 0.55],[0.26, 0, 0.73]])
    
    print("Test:")
    print("-----")
    
    print("[b,B]=AcyclicDPHFromMG(a,A):")
    b,B=AcyclicDPHFromMG(a,A)
    
    print("Moments of (a,A) and (b,B)")
    ma=MomentsFromMG(a,A,5)
    mb=MomentsFromMG(b,B,5)
    print(ma)
    print(mb)

    print("Check the obtained DPH, CheckDPHRepresentation(b,B):")
    flag=CheckDPHRepresentation(b,B)
    print(flag)     
   
    C=SimilarityMatrix(A,B)
    print(A*C)
    err1 = np.max(np.max(np.abs(A*C-C*B)))
    err2 = np.max(np.abs(a*C-b))
    
    print("Transformation errors:")
    print(err1)
    print(err2)
    
    assert err1<1e-12 and err2<1e-12, "Transformation to canonical PH(3) failed!"
       
    print("--DPH2From3Moments-------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.9, 0.1]])
    A=ml.matrix([[0.2, 0.61],[0.58, 0.41]])
    moms=MomentsFromDPH(a,A)
    print(moms)        
    
    print("Test:")
    print("-----")
    
    print("[a,A]=DPH2From3Moments(moms)")
    b,B=DPH2From3Moments(moms)
    
    print("MomentsFromMG(a,A,3):")
    phmoms = MomentsFromMG(b,B,3)
    print(phmoms)
    
    assert np.all(np.abs((np.array(phmoms)-np.array(moms))/np.array(moms))<1e-12), "The moments of the result of DPH2From3Moments do not match the input!"
    
    print("--DPH3From5Moments-------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.7, 0.1, 0.2]])
    A=ml.matrix([[0.2, 0.51, 0.1],[0.58, 0.41, 0],[0.1, 0.4, 0.3]])
    moms=MomentsFromDPH(a,A)
    print(moms)        
    
    print("Test:")
    print("-----")
    
    print("[a,A]=DPH3From5Moments(moms)")
    b,B=DPH3From5Moments(moms)
    
    print("MomentsFromMG(a,A,5):")
    phmoms = MomentsFromMG(b,B,5)
    print(phmoms)
    
    assert np.all(np.abs((np.array(phmoms)-np.array(moms))/np.array(moms))<1e-12), "The moments of the result of DPH3From5Moments do not match the input!"
    
    print("--SamplesFromDPH---------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.76, 0, 0.24]])
    print(a)
    A=ml.matrix([[0.34, 0.66, 0],[0.79, 0.05, 0.07],[0.26, 0.73, 0.01]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("x=SamplesFromDPH(a,A,1000)")
    x=SamplesFromDPH(a,A,1000)
    
    print("Moments from the samples:")
    mt = MarginalMomentsFromTrace(x,3)
    print(mt)
    
    print("Moments from the DPH:")
    mp = MomentsFromDPH(a,A,3)
    print(mp)

if __name__ == "__main__":
    TestDPHPackage()