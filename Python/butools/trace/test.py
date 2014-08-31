# -*- coding: utf-8 -*-
"""
Created on Sun Aug 31 15:45:54 2014

@author: gabor
"""

import butools
import numpy as np
import numpy.matlib as ml
from butools.moments import *

def TestTracePackage ():

    print("---BuTools: Trace package test file---")

    print("Enable the verbose messages with the BuToolsVerbose flag")
    butools.verbose = True
    
    print("Enable input parameter checking with the BuToolsCheckInput flag")
    butools.checkInput = True
   
    print("----------------------------------------------------------------------------")
    print("Generating trace file for the tests...")
    
    D0 = ml.matrix([[-18, 1, 4],[2, -18, 7],[1, 3, -32]])
    D1 = ml.matrix([[12, 1, 0],[1, 8, 0],[2, 1, 25]])
    tr = SamplesFromMAP(D0,D1,1000000)
    mmoms = MarginalMomentsFromMAP(D0,D1,5)
    macf = LagCorrelationsFromMAP(D0,D1,10)
    mNm1 = LagkJointMomentsFromMAP(D0,D1,3,1)
    mNm2 = LagkJointMomentsFromMAP(D0,D1,3,2)
    [a,A]=MarginalDistributionFromMAP(D0,D1)
    
    print("--CdfFromTrace----------------------------------------------------------------------")
    
    print("Test:")
    print("-----")
    
    print("[x,y]=CdfFromTrace(tr):")
    [x,y]=CdfFromTrace(tr)

    meandiff = np.abs(np.diff(x).dot(1.0-y[:-2]) - np.mean(tr))/np.mean(tr)
    assert np.all(np.diff(y)>=0) and y[-1]<=1 and y[0]>=0, "CdfFromTrace returned a wrong cdf!"
    assert meandiff<1e-2, "The mean obtained from the cdf returned by CdfFromTrace does not match the trace mean!"
    
    print("--PdfFromTrace----------------------------------------------------------------------")
    
    print("Test:");
    print("-----");
    
    print("[x,y]=PdfFromTrace(tr, 0:0.01:0.5):")
    [x,y]=PdfFromTrace(tr, 0:0.01:0.5)
    [xm,ym]=IntervalPdfFromPH(a, A, 0:0.01:0.5)
    
    meandiff = np.abs(x[0:-2].dot(np.diff(x)*y[0:-2]) - np.mean(tr))/np.mean(tr)
    
    assert np.all(y>=0), "PdfFromTrace returned a wrong pdf!"
    assert meandiff<1e-2, "The mean obtained from the pdf returned by PdfFromTrace does not match the trace mean!"
    
    print("----------------------------------------------------------------------------");
    help MarginalMomentsFromTrace
    
    print("Test:");
    print("-----");
    
    print("moms=MarginalMomentsFromTrace(tr, 3):");
    moms=MarginalMomentsFromTrace(tr, 3);
    moms
    mmoms(1:3)
    
    momdiff = norm((mmoms(1:3)-moms)./mmoms(1:3));
    
    assert(momdiff<1e-1, "Moments from MarginalMomentsFromTrace are far from the theoretical moments of the trace!");
    
    print("----------------------------------------------------------------------------");
    help LagCorrelationsFromTrace
    
    print("Test:");
    print("-----");
    
    print("acf=LagCorrelationsFromTrace(tr, 10):");
    acf=LagCorrelationsFromTrace(tr, 10);
    acf
    macf
    plot([acf,macf]);
    
    acfdiff = norm(acf-macf);
    
    assert(acfdiff<1e-1, "Autocorrelations from LagCorrelationsFromTrace are far from the theoretical autocorrelations of the trace!");
    
    print("----------------------------------------------------------------------------");
    help LagkJointMomentsFromTrace
    
    
    print("Test:");
    print("-----");
    
    print("Nm=LagkJointMomentsFromTrace(tr,3,1):");
    Nm=LagkJointMomentsFromTrace(tr,3,1)
    mNm1
    
    Nmdiff = norm((Nm-mNm1)./mNm1);
    assert(Nmdiff<1e-1, "Moments from LagkJointMomentsFromTrace are far from the theoretical moments of the trace!");
    
    print("Test:");
    print("-----");
    
    print("Nm=LagkJointMomentsFromTrace(tr,3,2):");
    Nm=LagkJointMomentsFromTrace(tr,3,2)
    mNm2
    
    Nmdiff = norm((Nm-mNm2)./mNm2);
    assert(Nmdiff<1e-1, "Moments from LagkJointMomentsFromTrace are far from the theoretical moments of the trace!");
    
    print("----------------------------------------------------------------------------");
    help CdfFromWeightedTrace
    
    print("Input:");
    print("------");
    wtr=[0.12; 1.23; 0.546; 0.6765; 1.34; 2.34];
    wei=[12; 1; 34; 23; 8; 2];
    
    print("Test:");
    print("-----");
    
    print("[x,y]=CdfFromWeightedTrace(wtr,wei):");
    [x,y]=CdfFromWeightedTrace(wtr,wei);
    plot(x,y);
    
    assert(all(diff(y)>=0) && y(end)<=1 && y(1)>=0, "CdfFromWeightedTrace returned a wrong cdf!");
    
    print("----------------------------------------------------------------------------");
    help PdfFromWeightedTrace
    
    print("Test:");
    print("-----");
    
    print("[x,y]=PdfFromWeightedTrace(wtr, wei 0:0.1:3):");
    [x,y]=PdfFromWeightedTrace(wtr, wei, 0:0.1:3);
    plot(x,y);
    
    assert(all(y>=0), "PdfFromWeightedTrace returned a wrong pdf!");
    
    print("----------------------------------------------------------------------------");
    help MarginalMomentsFromWeightedTrace
    
    print("Test:");
    print("-----");
    
    print("moms=MarginalMomentsFromWeightedTrace(wtr, wei, 3):");
    moms=MarginalMomentsFromWeightedTrace(wtr, wei, 3);
    moms
    

if __name__ == "__main__":
    TestTracePackage()