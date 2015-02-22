# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 17:20:24 2014

@author: gabor
"""

import numpy as np
import butools
from butools.map import *
from butools.trace import *
from butools.fitting import *
from butools.ph import *
import zipfile

def TestFittingPackage ():

    print("---BuTools: Fitting package test file---")
    
    print("Enable the verbose messages with the BuToolsVerbose flag")
    butools.verbose = True
    
    print("Enable input parameter checking with the BuToolsCheckInput flag")
    butools.checkInput = True
    

    zip = zipfile.ZipFile('../Matlab/fitting/bctrace.zip')  
    zip.extractall()  
    tr = np.loadtxt('bctrace.iat')
    
    print('Length of the trace:',len(tr),'samples')

    print('----------------------------------------------------------------------------')
    
    print("Generating input for the test")
    print("=============================")
    
    intBounds = np.linspace(0, MarginalMomentsFromTrace(tr,1)[0]*4, 50)
    print("Obtaining pdf of the trace: [pdfTrX, pdfTrY] = PdfFromTrace(tr,intBounds)")
    pdfTrX, pdfTrY = PdfFromTrace(tr,intBounds)
    print("Obtaining cdf of the trace: [cdfTrX, cdfTrY] = CdfFromTrace(tr)")
    cdfTrX, cdfTrY = CdfFromTrace(tr)
    
    maxCdfPoints = 2000
    if len(tr)>maxCdfPoints:
        step = math.ceil(len(tr) / maxCdfPoints)
        cdfTrX = cdfTrX[0:len(tr):step]
        cdfTrY = cdfTrY[0:len(tr):step]
    
    print("Fitting an APH based on 3 moments to the trace...")
    moms = MarginalMomentsFromTrace(tr,3)
    [alpha,A] = APHFrom3Moments(moms)
    print("Obtaining pdf of the APH(3): [pdfPHX, pdfPHY] = IntervalPdfFromPH(alpha, A, intBounds)")
    pdfPHX, pdfPHY = IntervalPdfFromPH(alpha, A, intBounds)
    print("Obtaining cdf of the APH(3): cdfPHY = CdfFromPH(alpha, A, cdfTrX)")
    cdfPHY = CdfFromPH(alpha, A, cdfTrX)
    
    print("----------------------------------------------------------------------------")

    print("Testing distance functions")
    print("==========================")
    
    print("Calculating pdf squared difference: sqPdf = EmpiricalSquaredDifference (pdfTrY, pdfPHY, intBounds)")
    sqPdf = EmpiricalSquaredDifference (pdfTrY, pdfPHY, intBounds)
    print(sqPdf)
    
    print("Calculating pdf relative entropy: rePdf = EmpiricalRelativeEntropy (pdfTrY, pdfPHY, intBounds)")
    rePdf = EmpiricalRelativeEntropy (pdfTrY, pdfPHY, intBounds)
    print(rePdf)
    
    print("Calculating cdf squared difference: sqCdf = EmpiricalSquaredDifference (cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX)")
    sqCdf = EmpiricalSquaredDifference (cdfTrY[0:-1], cdfPHY[0:-1], cdfTrX)
    print(sqCdf)
    
    print("Calculating cdf relative entropy: reCdf = EmpiricalRelativeEntropy (cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX)")
    reCdf = EmpiricalRelativeEntropy (cdfTrY[0:-1], cdfPHY[0:-1], cdfTrX)
    print(reCdf)
    
    print("Calculating likelihood: logli = LikelihoodFromTrace(tr,alpha,A)")   
    logli = LikelihoodFromTrace(tr,alpha,A)
    print(logli)
    
    print("----------------------------------------------------------------------------")

    print("Fitting a PH(5) by using G-FIT")
    print("==============================")
    
    print("Perform fitting: [alpha,A]=PHFromTrace(tr,5)")   
    [alpha,A]=PHFromTrace(tr,5)
    
    pdfPHX, pdfPHY = IntervalPdfFromPH(alpha, A, intBounds)
    cdfPHY = CdfFromPH(alpha, A, cdfTrX)

    print("Calculating pdf squared difference: sqPdf = EmpiricalSquaredDifference (pdfTrY, pdfPHY, intBounds)")
    sqPdf = EmpiricalSquaredDifference (pdfTrY, pdfPHY, intBounds)
    print(sqPdf)
    
    print("Calculating pdf relative entropy: rePdf = EmpiricalRelativeEntropy (pdfTrY, pdfPHY, intBounds)")
    rePdf = EmpiricalRelativeEntropy (pdfTrY, pdfPHY, intBounds)
    print(rePdf)
    
    print("Calculating cdf squared difference: sqCdf = EmpiricalSquaredDifference (cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX)")
    sqCdf = EmpiricalSquaredDifference (cdfTrY[0:-1], cdfPHY[0:-1], cdfTrX)
    print(sqCdf)
    
    print("Calculating cdf relative entropy: reCdf = EmpiricalRelativeEntropy (cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX)")
    reCdf = EmpiricalRelativeEntropy (cdfTrY[0:-1], cdfPHY[0:-1], cdfTrX)
    print(reCdf)
    
    print("Calculating likelihood: logli = LikelihoodFromTrace(tr,alpha,A)")   
    logli = LikelihoodFromTrace(tr,alpha,A)
    print(logli)
    
    print("----------------------------------------------------------------------------")

    print("Testing further distance functions")
    print("==================================")

    maxMAPTraceLen = 10000
    if len(tr)>maxMAPTraceLen:
        tr = tr[0:maxMAPTraceLen]

    print("Fitting a MAP based on 3 moments and 1 lag correlation...")
    corr1 = LagCorrelationsFromTrace(tr,1)[0]
    [D0,D1]=MAPFromFewMomentsAndCorrelations(moms,corr1)

    print("Calculating likelihood: logli = LikelihoodFromTrace(tr,D0,D1)")   
    logli = LikelihoodFromTrace(tr,D0,D1)
    print(logli)
    
    print("Obtaining acf of the trace: trAcf = LagCorrelationsFromTrace(tr, 10)")
    trAcf = LagCorrelationsFromTrace(tr, 10)
    print(trAcf)
    
    print("Obtaining acf of the MAP: mAcf = LagCorrelationsFromMAP(D0, D1, 10)")
    mAcf = LagCorrelationsFromMAP(D0, D1, 10, 1e-13)
    print(mAcf)

    print("Calculating acf squared difference: sqAcf = SquaredDifference (mAcf, trAcf)")
    sqAcf = SquaredDifference (mAcf, trAcf)
    print(sqAcf)
    
    print("Calculating acf relative entropy: reAcf = RelativeEntropy (mAcf, trAcf)")
    reAcf = RelativeEntropy (mAcf, trAcf)
    print(reAcf)
    
    print("Fitting a MAP(5) by using SPEM-FIT")
    print("==================================")
    
    print("Perform fitting: [D0,D1]=MAPFromTrace(tr,5)")
    [D0,D1]=MAPFromTrace(tr,5)
    
    print("Calculating likelihood: logli = LikelihoodFromTrace(tr,D0,D1)")   
    logli = LikelihoodFromTrace(tr,D0,D1)
    print(logli)

    print("Obtaining acf of the MAP: mAcf = LagCorrelationsFromMAP(D0, D1, 10)")
    mAcf = LagCorrelationsFromMAP(D0, D1, 10, 1e-13)
    print(mAcf)

    print("Calculating acf squared difference: sqAcf = SquaredDifference (mAcf, trAcf)")
    sqAcf = SquaredDifference (mAcf, trAcf)
    print(sqAcf)
    
    print("Calculating acf relative entropy: reAcf = RelativeEntropy (mAcf, trAcf)")
    reAcf = RelativeEntropy (mAcf, trAcf)
    print(reAcf)

if __name__ == "__main__":
    TestFittingPackage()