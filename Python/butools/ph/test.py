# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 17:20:24 2014

@author: gabor
"""

import numpy as np
import numpy.matlib as ml
import butools
from butools.ph import *
from butools.moments import CheckMoments
from butools.trace import MarginalMomentsFromTrace

def TestPHPackage ():

    print('---BuTools: PH package test file---')
    
    print('Enable the verbose messages with the BuToolsVerbose flag')
    butools.verbose = True
    
    print('Enable input parameter checking with the BuToolsCheckInput flag')
    butools.checkInput = True
    
    print("--MomentsFromME---------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a = ml.matrix([[0.2, 0.3, 0.5]])
    print(a)
    A = ml.matrix([[-1,0,0],[0,-3,2],[0,-2,-3]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("MomentsFromME(a,A):")
    moms=MomentsFromME(a,A)
    print(moms)
    
    assert CheckMoments(moms)==True, "The function returned invalid moments!"
    
    print("Test:")
    print("-----")
    
    print("MomentsFromME(a,A,9):")
    moms=MomentsFromME(a,A,9)
    print(moms)
    
    assert CheckMoments(moms)==True, "The function returned invalid moments!"
    
    print("--MomentsFromPH---------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.1, 0.9, 0]])
    print(a)
    A=ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("MomentsFromPH(a,A,5)")
    moms=MomentsFromPH(a,A,5)
    print(moms)
    
    assert CheckMoments(moms)==True, "The function returned invalid moments!"
    
    print("--CdfFromME-------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a = ml.matrix([[0.2, 0.3, 0.5]])
    print(a)
    A = ml.matrix([[-1,0,0],[0,-3,2],[0,-2,-3]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("cdf=CdfFromME(a,A,(0:0.01:5)'):")
    cdf = CdfFromME(a,A,np.linspace(0,5,500))
    
    assert np.all(np.diff(cdf)>0), "The cdf is not increasing monotonously!"
    assert abs(np.sum(1.0-cdf)*0.01 - MomentsFromME(a,A,1))<0.01, "The mean computed from the cdf does not match the theoretical result!"
    
    print("--CdfFromPH-------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.1, 0.9, 0]])
    print(a)
    A=ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("cdf=CdfFromPH(a,A,(0:0.002:3)''):")
    cdf = CdfFromPH(a,A,np.linspace(0,3,1500))
    
    assert np.all(np.diff(cdf)>0), "The cdf is not increasing monotonously!"
    assert abs(np.sum(1.0-cdf)*0.002 - MomentsFromPH(a,A,1))<0.01, "The mean computed from the cdf does not match the theoretical result!"
    
    print("--PdfFromME-------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a = ml.matrix([[0.2, 0.3, 0.5]])
    print(a)
    A = ml.matrix([[-1,0,0],[0,-3,2],[0,-2,-3]])
    print(A)
    x = np.linspace(0,5,500);
    
    print("Test:")
    print("-----")
    
    print("pdf=PdfFromME(a,A,x):")
    pdf = PdfFromME(a,A,x)
    
    assert np.all(pdf)>=0, "The pdf is negative!"
    assert abs(np.dot(pdf,x)*0.01 - MomentsFromME(a,A,1))<0.01, "The mean computed from the pdf does not match the theoretical result!"
    
    print("--PdfFromPH-------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.1, 0.9, 0]])
    print(a)
    A=ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    print(A)
    x = np.linspace(0,3,1500)
    
    print("Test:")
    print("-----")
    
    print("pdf=PdfFromPH(a,A,x):")
    pdf = PdfFromPH(a,A,x)
    
    assert np.all(pdf)>=0, "The pdf is negative!"
    assert abs(np.dot(pdf,x)*0.002 - MomentsFromPH(a,A,1))<0.002, "The mean computed from the pdf does not match the theoretical result!"
    
    print("--IntervalPdfFromPH-----------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.1, 0.9, 0]])
    print(a)
    A=ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    print(A)
    x = np.linspace(0,3,1500)
    
    print("Test:")
    print("-----")
    
    print("[x,y]=IntervalPdfFromPH(a,A,x):")
    [x,y] = IntervalPdfFromPH(a,A,x)
       
    assert np.all(y)>=0, "The interval pdf is negative!"
    assert abs(np.dot(y,x)*0.002 - MomentsFromPH(a,A,1))<0.01, "The mean computed from the interval pdf does not match the theoretical result!"
    
    print("--RandomPH--------------------------------------------------------------------------")
    
    print("Test:")
    print("-----")
       
    print("[a,A]=RandomPH(3,8,4):")
    a,A=RandomPH(3,8,4)
       
    assert CheckPHRepresentation(a,A), "RandomPH failed to return a valid PH representation!"
    assert np.max(np.abs(MomentsFromPH(a,A,1)[0]-8.0))<1e-14, "RandomPH failed to match the given mean value!"
    
    print("--CheckMERepresentation-------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[-0.2, 0.2]])
    print(a)
    A=ml.matrix([[1, -1],[1, -2]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("CheckMERepresentation(a,A):")
    flag=CheckMERepresentation(a,A)
    print(flag)
    
    assert flag==False, "CheckMERepresentation did not detect that the initial vector is invalid!"
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[-0.2, 0.4, 0.8]])
    print(a)
    A=ml.matrix([[-2, 0, 3],[0, -1, 1],[0, -1, -1]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("CheckMERepresentation(a,A):")
    flag=CheckMERepresentation(a,A)
    print(flag)
    
    assert flag==False, "CheckMERepresentation did not detect that the dominant eigenvalue is invalid!"
    
    print("Input:")
    print("------")
    
    a = ml.matrix([[0.2, 0.3, 0.5]])
    print(a)
    A = ml.matrix([[-1,0,0],[0,-3,2],[0,-2,-3]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("CheckMERepresentation(a,A):")
    flag=CheckMERepresentation(a,A)
    print(flag)
    
    assert flag==True, "CheckMERepresentation did not recognize that the given ME representation is valid!"
    
    print("--CheckPHRepresentation-------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.2]])
    A=ml.matrix([[-1, 1],[1,-2]])
    
    print("Test:")
    print("-----")
    
    print("CheckPHRepresentation(a,A):")
    flag=CheckPHRepresentation(a,A)
    print(flag)
    
    assert flag==False, "CheckPHRepresentation did not recognize the wrong input dimensions!"
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.2, 0.7]])
    A=ml.matrix([[-1, 1],[1, -2]])
    
    print("Test:")
    print("-----")
    
    print("CheckPHRepresentation(a,A):")
    flag=CheckPHRepresentation(a,A)
    print(flag)
    
    assert flag==True, "CheckPHRepresentation did not recognize that the given PH representation is valid!"
    
    print("--CheckMEPositiveDensity------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a = ml.matrix([[0.2, 0.3, 0.5]])
    print(a)
    A = ml.matrix([[-1,0,0],[0,-3,2],[0,-2,-3]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("CheckMEPositiveDensity(a,A):")
    flag=CheckMEPositiveDensity(a,A)
    print(flag)
    
    assert flag==True, "CheckMEPositiveDensity did not recognize that the given ME distribution has positive density!"
    
    print("Input:")
    print("------")
    
    a = ml.matrix([[0.2, 0.3, 0.5]])
    print(a)
    A = ml.matrix([[-1,0,0],[0,-3,2.9],[0,-2.9,-3]])
    print(A)

    print("Test:")
    print("-----")
    
    print("CheckMEPositiveDensity(a,A):")
    flag=CheckMEPositiveDensity(a,A)
    print(flag)
    
    assert flag==False, "CheckMEPositiveDensity did not recognize that the given ME distribution does not have positive density!"
    
    print("--APHFrom3Moments-------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    moms = [10,125,8400]
    
    print("Test:")
    print("-----")
    
    print("[a,A]=APHFrom3Moments(moms(1), moms(2), moms(3))")
    print("size(A,1)")
    a,A=APHFrom3Moments(moms)
    print(A.shape[0])
    
    print("MomentsFromPH(a,A,3):")
    phmoms = MomentsFromPH(a,A,3)
    print(phmoms)
    
    assert np.all(np.abs((np.array(phmoms)-np.array(moms))/np.array(moms))<1e-12), "APHFrom3Moments failed to match the given moments!"
    
    print("Input:")
    print("------")
    
    moms = [10,525,31400]
    
    print("Test:")
    print("-----")
    
    print("[a,A]=APHFrom3Moments(moms(1), moms(2), moms(3))")
    print("size(A,1)")
    a,A=APHFrom3Moments(moms)
    print(A.shape[0])
    
    print("MomentsFromPH(a,A,3):")
    phmoms = MomentsFromPH(a,A,3)
    print(phmoms)
    
    assert np.all(np.abs((np.array(phmoms)-np.array(moms))/np.array(moms))<1e-12), "APHFrom3Moments failed to match the given moments!"
    
    
    print("--PH2From3Moments-------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    moms = [10,160,3500]
    
    print("Test:")
    print("-----")
    
    print("[a,A]=PH2From3Moments(moms)")
    a,A=PH2From3Moments(moms)
    
    print("MomentsFromPH(a,A,3):")
    phmoms = MomentsFromPH(a,A,3)
    print(phmoms)
    
    assert np.all(np.abs((np.array(phmoms)-np.array(moms))/np.array(moms))<1e-12), "PH2From3Moments failed to match the given moments!"
    
    print("Input:")
    print("------")
    
    moms = [10,260,13500]
    
    print("Test:")
    print("-----")
    
    print("[a,A]=PH2From3Moments(moms)")
    a,A=PH2From3Moments(moms)
    
    print("MomentsFromPH(a,A,3):")
    phmoms = MomentsFromPH(a,A,3)
    print(phmoms)
    
    assert np.all(np.abs((np.array(phmoms)-np.array(moms))/np.array(moms))<1e-12), "PH2From3Moments failed to match the given moments!"
    
    print("--PH3From5Moments-------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.1, 0.9, 0]])
    print(a)
    A=ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    print(A)

    moms = MomentsFromPH(a,A)
    print(moms)
    
    print("Test:")
    print("-----")
    
    print("[a,A]=PH3From5Moments(moms)")
    a,A=PH3From5Moments(moms)
    
    print("MomentsFromPH(a,A,5):")
    phmoms = MomentsFromPH(a,A,5)
    print(phmoms)
    
    assert np.all(np.abs((np.array(phmoms)-np.array(moms))/np.array(moms))<1e-12), "PH3From5Moments failed to match the given moments!"
    
    print("Input:")
    print("------")
    
    a = ml.matrix([[0.2, 0.3, 0.5]])
    print(a)
    A = ml.matrix([[-1,0,0],[0,-3,0.5],[0,-0.5,-3]])
    print(A)
    moms = MomentsFromME(a,A)
    print(moms)
    
    print("Test:")
    print("-----")
    
    print("[a,A]=PH3From5Moments(moms)")
    a,A=PH3From5Moments(moms)
    
    print("MomentsFromPH(a,A,5):")
    phmoms = MomentsFromPH(a,A,5)
    print(phmoms)
    
    assert np.all(np.abs((np.array(phmoms)-np.array(moms))/np.array(moms))<1e-12), "PH3From5Moments failed to match the given moments!"
    
    print("--MEFromMoments---------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.1, 0.9, 0]])
    print(a)
    A=ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    print(A)

    moms=MomentsFromPH(a,A,5)
    print(moms)
    
    print("Test:")
    print("-----")
    
    print("[a,A]=MEFromMoments(moms):")
    a,A=MEFromMoments(moms)
    
    print("MomentsFromME(a,A,5):")
    memoms = MomentsFromME(a,A,5)
    print(memoms)
    
    assert np.all(np.abs((np.array(memoms)-np.array(moms))/np.array(moms))<1e-12), "MEFromMoments failed to match the given moments!"
    
    print("--APH2ndMomentLowerBound------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    mean = 1.9
    n = 4
    
    print("Test:")
    print("-----")
    
    print("mom2 = APH2ndMomentLowerBound(mean,n):")
    mom2 = APH2ndMomentLowerBound(mean,n)
    print(mom2)
    
    cv2 = mom2/mean**2-1
    
    assert abs(cv2-1/n)<1e-14, "APH2ndMomentLowerBound did not give the expected result!"
    
    
    print("--APH3rdMomentLowerBound/APH3rdMomentUpperBound-------------------------------------")
    
    print("Input:")
    print("------")
    
    mean = 1.9
    mom2 = 5
    n = 3
    
    print("Test:")
    print("-----")
    
    print("mom3lower = APH3rdMomentLowerBound(mean,mom2,n):")
    mom3lower = APH3rdMomentLowerBound(mean,mom2,n)
    print(mom3lower)
    print("mom3upper = APH3rdMomentUpperBound(mean,mom2,n):")
    mom3upper = APH3rdMomentUpperBound(mean,mom2,n)
    print(mom3upper)
    
    assert mom3upper>mom3lower, "Lower bound is larger than the upper bound!"
    
    print("Input:")
    print("------")
    
    mean = 1.9
    mom2 = 5
    n = 4
    
    print("Test:")
    print("-----")
    
    print("mom3lower = APH3rdMomentLowerBound(mean,mom2,n):")
    mom3lower = APH3rdMomentLowerBound(mean,mom2,n)
    print(mom3lower)
    print("mom3upper = APH3rdMomentUpperBound(mean,mom2,n):")
    mom3upper = APH3rdMomentUpperBound(mean,mom2,n)
    print(mom3upper)
    
    assert mom3upper>mom3lower, "Lower bound is larger than the upper bound!"
    assert mom3upper==np.inf, "Upper bound must be infinity with 4 phases!"
    
    print("--CanonicalFromPH2------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.12, 0.88]])
    A=ml.matrix([[-1.28, 0],[3.94, -3.94]])
    
    print("Test:")
    print("-----")
    
    print("[b,B]=CanonicalFromPH2(a,A):")
    b,B=CanonicalFromPH2(a,A)
    
    C=SimilarityMatrix(A,B)
    err1 = np.max(np.max(np.abs(A*C-C*B)))
    err2 = np.max(np.abs(a*C-b))
    
    print("Transformation errors:")
    print(err1)
    print(err2)
    
    assert err1<1e-12 and err2<1e-12, "Transformation to canonical PH(2) failed!"
    
    print("--CanonicalFromPH3------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.1, 0.9, 0]])
    print(a)
    A=ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("[b,B]=CanonicalFromPH3(a,A):")
    b,B=CanonicalFromPH3(a,A)
   
    C=SimilarityMatrix(A,B)
    print(A*C)
    err1 = np.max(np.max(np.abs(A*C-C*B)))
    err2 = np.max(np.abs(a*C-b))
    
    print("Transformation errors:")
    print(err1)
    print(err2)
    
    assert err1<1e-12 and err2<1e-12, "Transformation to canonical PH(3) failed!"
    
    print("--AcyclicPHFromME-------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[-0.4, 1.4, 0]])
    A=ml.matrix([[-4, 1, 1],[0, -2, 1],[1, 0, -8]])
    
    print("Test:")
    print("-----")
    
    print("[b,B]=AcyclicPHFromME(a,A):")
    b,B=AcyclicPHFromME(a,A)
    
    print("Moments of (a,A) and (b,B)")
    ma=MomentsFromME(a,A,5)
    mb=MomentsFromME(b,B,5)
    print(ma)
    print(mb)
    
    assert la.norm((np.array(ma)-np.array(mb))/np.array(ma))<1e-7, "Transformation to acyclic representation failed!"
    
    print("--MonocyclicPHFromME----------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a = ml.matrix([[0.2, 0.3, 0.5]])
    print(a)
    A = ml.matrix([[-1,0,0],[0,-3,2],[0,-2,-3]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("[b,B]=MonocyclicPHFromME(a,A):")
    b,B=MonocyclicPHFromME(a,A)
    
    print("Moments of (a,A) and (b,B)")
    ma=MomentsFromME(a,A,5)
    mb=MomentsFromME(b,B,5)
    print(ma)
    print(mb)
    
    assert la.norm((np.array(ma)-np.array(mb))/np.array(ma))<1e-7, "Transformation to monocyclic representation failed!"
    
    print("--PHFromME--------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[-0.4, 1.4]])
    A=ml.matrix([[-3.8, 2],[2, -9]])
    
    print("Test:")
    print("-----")
    
    print("CheckMERepresentation(a,A):")
    flag=CheckMERepresentation(a,A)
    print(flag)
    
    print("CheckPHRepresentation(a,A):")
    flag=CheckPHRepresentation(a,A)
    print(flag)
    
    print("[b,B]=PHFromME(a,A)")
    b,B=PHFromME(a,A)
    print(b)
    print(B)
    
    print("Check the obtained PH, CheckPHRepresentation(b,B):")
    flag=CheckPHRepresentation(b,B)
    print(flag)
    
    C=SimilarityMatrix(A,B)
    err1 = np.max(np.max(np.abs(A*C-C*B)))
    err2 = np.max(np.abs(a*C-b))
    print("Transformation errors:")
    print(err1)
    print(err2)
    
    assert flag and err1<1e-12 and err2<1e-12, "Transformation to PH failed!"
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[-0.5, 1.5]])
    A=ml.matrix([[-3.8, 2],[2, -9]])
    
    print("Test:")
    print("-----")
    
    print("CheckMERepresentation(a,A):")
    flag=CheckMERepresentation(a,A)
    print(flag)
    
    print("CheckPHRepresentation(a,A):")
    flag=CheckPHRepresentation(a,A)
    print(flag)
    
    print("[b,B]=PHFromME(a,A)")
    [b,B]=PHFromME(a,A)
    print(b)
    print(B)
    
    print("Check the obtained PH, CheckPHRepresentation(b,B):")
    flag=CheckPHRepresentation(b,B)
    print(flag)
    
    C=SimilarityMatrix(A,B)
    err1 = np.max(np.max(np.abs(A*C-C*B)))
    err2 = np.max(np.abs(a*C-b))
    print("Transformation errors:")
    print(err1)
    print(err2)
    
    assert flag and err1<1e-12 and err2<1e-12, "Transformation to PH failed!"
    
    print("--MEOrder---------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]])/6.0
    A=ml.matrix([[-1, 0, 0, 0, 0, 0],[0.5, -2, 1, 0, 0, 0],[1, 0, -3, 1, 0, 0],[1, 0, 1, -4, 1, 0],[4, 0, 0, 0, -5, 0],[5, 0, 0, 0, 0, -6]])
    print(a)
    print(A)
    
    print("Test:")
    print("-----")
    
    print("co=MEOrder(a,A,''cont'')")
    co=MEOrder(a,A,"cont")
    print(co)
    print("oo=MEOrder(a,A,''obs'')")
    oo=MEOrder(a,A,"obs")
    print(oo)
    print("coo=MEOrder(a,A,''obscont'')")
    coo=MEOrder(a,A,"obscont")
    print(coo)
    print("coo=MEOrder(a,A,''moment'')")
    mo=MEOrder(a,A,"moment")
    print(mo)
    
    assert co==2, "Wrong controllability order returned!"
    assert oo==6, "Wrong observability order returned!"
    assert coo==2, "The minimum of the controllability and observability order is wrong!"
    assert mo==2, "Wrong moment order returned!"
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[2, 1]])/3
    A=ml.matrix([[-1, 1],[0, -3]])
    print(a)
    print(A)    
    
    print("Test:")
    print("-----")
    
    print("co=MEOrder(a,A,''cont'')")
    co=MEOrder(a,A,"cont")
    print(co)
    print("oo=MEOrder(a,A,''obs'')")
    oo=MEOrder(a,A,"obs")
    print(oo)
    print("coo=MEOrder(a,A,''obscont'')")
    coo=MEOrder(a,A,"obscont")
    print(coo)
    print("coo=MEOrder(a,A,''moment'')")
    mo=MEOrder(a,A,"moment")
    print(mo)
    
    assert co==2, "Wrong controllability order returned!"
    assert oo==1, "Wrong observability order returned!"
    assert coo==1, "The minimum of the controllability and observability order is wrong!"
    assert mo==1, "Wrong moment order returned!"
    
    print("Input:")
    print("------")
    
    b = ml.matrix([[0.2, 0.3, 0.5]])
    B = ml.matrix([[-1,0,0],[0,-3,1],[0,-1,-3]])
    [a,A] = MonocyclicPHFromME(b,B)
    print(a)
    print(A)
    
    print("Test:")
    print("-----")
    
    print("co=MEOrder(a,A,''cont'')")
    co=MEOrder(a,A,"cont")
    print(co)
    print("oo=MEOrder(a,A,''obs'')")
    oo=MEOrder(a,A,"obs")
    print(oo)
    print("coo=MEOrder(a,A,''obscont'')")
    coo=MEOrder(a,A,"obscont")
    print(coo)
    print("coo=MEOrder(a,A,''moment'')")
    mo=MEOrder(a,A,"moment")
    print(mo)
    
    assert co==9, "Wrong controllability order returned!"
    assert oo==3, "Wrong observability order returned!"
    assert coo==3, "The minimum of the controllability and observability order is wrong!"
    assert mo==3, "Wrong moment order returned!"
    
    print("--MEOrderFromMoments----------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.1, 0.9, 0]])
    print(a)
    A=ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    print(A)
    
    print("Test:")
    print("-----")
    
    print("moms=MomentsFromME(a,A)")
    print("mo=MEOrderFromMoments(moms):")
    moms=MomentsFromME(a,A)
    print(moms)
    mo = MEOrderFromMoments(moms)
    print(mo)
    
    assert mo==3, "Wrong moment order returned!"
    
    print("Input:")
    print("------")
    
    b = ml.matrix([[0.2, 0.3, 0.5]])
    B = ml.matrix([[-1,0,0],[0,-3,2],[0,-2,-3]])
    [a,A] = MonocyclicPHFromME(b,B)
    print(a)
    print(A)
    
    print("Test:")
    print("-----")
    
    print("moms=MomentsFromME(a,A)")
    print("mo=MEOrderFromMoments(moms):")
    moms=MomentsFromME(a,A)
    print(moms)
    mo = MEOrderFromMoments(moms)
    print(mo)
    
    assert mo==3, "Wrong moment order returned!"
    
    print("--MinimalRepFromME------------------------------------------------------------------")
    help 
    
    print("Input:")
    print("------")

    a=ml.matrix([[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]])/6.0
    A=ml.matrix([[-1, 0, 0, 0, 0, 0],[0.5, -2, 1, 0, 0, 0],[1, 0, -3, 1, 0, 0],[1, 0, 1, -4, 1, 0],[4, 0, 0, 0, -5, 0],[5, 0, 0, 0, 0, -6]])
    print(a)
    print(A)
    
    print("Test:")
    print("-----")
    
    print("[b,B]=MinimalRepFromME(a,A,''cont'')")
    [b,B]=MinimalRepFromME(a,A,"cont")
    print(b)
    print(B)
    
    assert b.shape[1]==2, "Non-minimal representation returned based on controllability!"
    
    print("[b,B]=MinimalRepFromME(a,A,''obs'')")
    [b,B]=MinimalRepFromME(a,A,"obs")
    print(b)
    print(B)
    
    assert b.shape[1]==6, "Non-minimal representation returned based on observability!"
    
    print("[b,B]=MinimalRepFromME(a,A,''obscont'')")
    [b,B]=MinimalRepFromME(a,A,"obscont")
    print(b)
    print(B)
    
    assert b.shape[1]==2, "Non-minimal representation returned based on observability and controllability!"
    
    print("[b,B]=MinimalRepFromME(a,A,''moment'')")
    [b,B]=MinimalRepFromME(a,A,"moment")
    print(b)
    print(B)
    
    assert b.shape[1]==2, "Non-minimal representation returned based on the moments!"
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[2, 1]])/3
    A=ml.matrix([[-1, 1],[0, -3]])
    print(a)
    print(A)
    
    print("Test:")
    print("-----")
    
    print("[b,B]=MinimalRepFromME(a,A,''cont'')")
    [b,B]=MinimalRepFromME(a,A,"cont")
    print(b)
    print(B)
    
    assert b.shape[1]==2, "Non-minimal representation returned based on controllability!"
    
    print("[b,B]=MinimalRepFromME(a,A,''obs'')")
    [b,B]=MinimalRepFromME(a,A,"obs")
    print(b)
    print(B)
    
    assert b.shape[1]==1, "Non-minimal representation returned based on observability!"
    
    print("[b,B]=MinimalRepFromME(a,A,''obscont'')")
    [b,B]=MinimalRepFromME(a,A,"obscont")
    print(b)
    print(B)
    
    assert b.shape[1]==1, "Non-minimal representation returned based on observability and controllability!"
    
    print("[b,B]=MinimalRepFromME(a,A,''moment'')")
    [b,B]=MinimalRepFromME(a,A,"moment")
    print(b)
    print(B)
    
    assert b.shape[1]==1, "Non-minimal representation returned based on the moments!"
    
    print("Input:")
    print("------")
    
    b = ml.matrix([[0.2, 0.3, 0.5]])
    B = ml.matrix([[-1,0,0],[0,-3,1],[0,-1,-3]])
    [a,A] = MonocyclicPHFromME(b,B)
    print(a)
    print(A)
    
    print("[b,B]=MinimalRepFromME(a,A,''cont'')")
    [b,B]=MinimalRepFromME(a,A,"cont")
    print(b)
    print(B)
    
    assert b.shape[1]==a.shape[1], "Non-minimal representation returned based on controllability!"
    
    print("[b,B]=MinimalRepFromME(a,A,''obs'')")
    [b,B]=MinimalRepFromME(a,A,"obs")
    print(b)
    print(B)
    # check similarity:
    C = SimilarityMatrix(B,A)
    sim = la.norm(B*C-C*A) + la.norm(b*C-a)
    print(sim)
    
    assert b.shape[1]==3 and sim<1e-13, "Non-minimal representation returned based on observability!"
    
    print("[b,B]=MinimalRepFromME(a,A,''obscont'')")
    [b,B]=MinimalRepFromME(a,A,"obscont")
    print(b)
    print(B)
    # check similarity:
    C = SimilarityMatrix(B,A)
    sim = la.norm(B*C-C*A) + la.norm(b*C-a)
    print(sim)
    
    assert b.shape[1]==3 and sim<1e-13, "Non-minimal representation returned based on observability and controllability!"
    
    print("[b,B]=MinimalRepFromME(a,A,''moment'')")
    [b,B]=MinimalRepFromME(a,A,"moment")
    print(b)
    print(B)
    # check similarity:
    C = SimilarityMatrix(B,A)
    sim = la.norm(B*C-C*A) + la.norm(b*C-a)
    print(sim)
    
    assert b.shape[1]==3 and sim<1e-13, "Non-minimal representation returned based on the moments!"
    
    print("--SamplesFromPH---------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    a=ml.matrix([[0.1, 0.9, 0]])
    A=ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    print(a)
    print(A)
    
    print("Test:")
    print("-----")
    
    print("x=SamplesFromPH(a,A,1000)")
    x=SamplesFromPH(a,A,1000)
    
    print("Moments from the samples:")
    mt = MarginalMomentsFromTrace(x,3)
    print(mt)
    
    print("Moments from the PH:")
    mp = MomentsFromPH(a,A,3)
    print(mp)

if __name__ == "__main__":
    TestPHPackage()