# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 17:20:24 2014

@author: gabor
"""

import sys
import numpy as np
import numpy.matlib as ml
import butools;
from butools.mc import *


def TestMCPackage ():

    print('---BuTools: MC package test file---')
    
    print('Enable the verbose messages with the BuToolsVerbose flag')
    butools.verbose = True
    
    print('Enable input parameter checking with the BuToolsCheckInput flag')
    butools.checkInput = True
    
    print('--CRPSolve--------------------------------------------------------------------------')
    
    print('Input:')
    print('------')
    
    Q1= ml.matrix("[-4.3 3.5 0.8; -8.4 6.5 1.9; 17.3 -12.7 -4.6]")
    print(Q1)
    
    print('Test:')
    print('-----')
    
    print('CRPSolve(Q1):')
    ret=CRPSolve(Q1)
    print(ret)
    
    assert np.abs(np.max(ret*Q1))<1e-14, 'The solution does not satisfy ret*Q1=0!'
    
    print('--DRPSolve--------------------------------------------------------------------------')
    
    print('Input:')
    print('------')
    
    Q2=ml.matrix("[-0.9 0.5 1.4; 0.9 -0.9 1; 0.3 1.3 -0.6]")
    print(Q2)
    
    print('Test:')
    print('-----')
    
    print('DRPSolve(Q2):')
    ret=DRPSolve(Q2)
    print(ret)
    assert np.abs(np.max(ret*Q2-ret))<1e-14, 'The solution does not satisfy ret*Q2=ret!'
    
    print('--CTMCSolve--------------------------------------------------------------------------')
    
    print('Input:')
    print('------')
    
    Q3 = ml.matrix("[0.1 0.5 0.4; 0.9 0.1 0]")
    print(Q3)
    
    print('Test:')
    print('-----')
    
    print('CTMCSolve(Q3):')
    try:
        ret=CTMCSolve(Q3)
        print(ret)
    except Exception as e:
        print(e)
    
    print('Input:')
    print('------')
    
    Q4 = ml.matrix("[0.1 0.5 0.4; 0.9 0.1 0; 0.9 0.1 0]")
    print(Q4)
    
    print('Test:')
    print('-----')
    
    print('CTMCSolve(Q4):')
    
    try:
        ret=CTMCSolve(Q4)
        print(ret)
    except Exception as e:
        print(e)
        
    print('Input:')
    print('------')
    
    Q5 = ml.matrix("[-0.9 0.5 0.4; 0.9 -0.9 0; 0.3 0.3 -0.6]")
    print(Q5)    
    
    print('Test:')
    print('-----')
    
    print('CTMCSolve(Q5):')
    ret=CTMCSolve(Q5)
    print(ret)
    assert np.abs(np.max(ret*Q5))<1e-14, 'The solution does not satisfy ret*Q5=0!'
    
    print('--DTMCSolve--------------------------------------------------------------------------')
    
    print('Input:')
    print('------')
    
    Q6 = ml.matrix("[0.1 0.5 0.4; 0.9 0.1 0; 0.3 -0.3 0.4]")
    print(Q6)    
    
    print('Test:')
    print('-----')
    
    print('DTMCSolve(Q6):')
    try:
        ret=DTMCSolve(Q6)
        print(ret)
    except Exception as e:
        print(e)
        
    print('Input:')
    print('------')
    
    Q7 = ml.matrix("[0.1 0.5 0.4; 0.9 0.1 0; 0.3 0.3 0.4]")
    print(Q7)
    
    print('Test:')
    print('-----')
    
    print('DTMCSolve(Q7):')
    ret=DTMCSolve(Q7)
    print(ret)
    assert np.abs(np.max(ret*Q7-ret))<1e-14, 'The solution does not satisfy ret*Q7=ret!'
    
    print('--CheckGenerator--------------------------------------------------------------------------')
    
    print('Input:')
    print('------')
    
    Q8 = ml.matrix("[-0.9 0.2 0.4; 0 0.9 0.9; 0 0.6 -0.6]")
    print(Q8)    
    
    print('Test:')
    print('-----')
    
    print('CheckGenerator(Q8,True):')
    flag=CheckGenerator(Q8,True)
    print(flag)
    assert flag==0,'CheckGenerator did not detect bad row sum!'
    
    print('Input:')
    print('------')
    
    Q9 = ml.matrix("[-0.9 0.5 0.4; 0.9 -0.9 0; 0.3 0.3 -0.6]")
    print(Q9)    
    
    print('Test:')
    print('-----')
    
    print('CheckGenerator(Q9,True):')
    flag=CheckGenerator(Q9,True)
    print(flag)
    assert flag==1,'CheckGenerator did not recognize a valid input!'
    
    
    print('Input:')
    print('------')
    
    Q10 = ml.matrix("[-0.9 0.2 0.4; 0.9 -0.9 0; 0.3 0.3 -0.6]")
    print(Q10)
    
    print('Test:')
    print('-----')
    
    print('CheckGenerator(Q10,True):')
    flag=CheckGenerator(Q10,True)
    print(flag)
    assert flag==1,'CheckGenerator did not recognize a valid input!'
    
    print('Input:')
    print('------')
    
    Q11 = ml.matrix("[-0.9 0.5 0.4; 0.9 -1.1 0; 0.3 0.3 -0.6]")
    print(Q11)    
    
    print('Test:')
    print('-----')
    
    print('CheckGenerator(Q11):')
    flag=CheckGenerator(Q11)
    print(flag)
    assert flag==0,'CheckGenerator did not recognize the non-zero row sum!'
    
    print('Input:')
    print('------')
    
    Q12 = ml.matrix("[-0.9 0.5 0.4; 0.9 -0.9 0; 0.3 0.3 -0.6]")
    print(Q12)    
    
    print('Test:')
    print('-----')
    
    print('CheckGenerator(Q12):')
    flag=CheckGenerator(Q12)
    print(flag)
    assert flag==1,'CheckGenerator did not recognize a valid input!'
    
    print('--CheckProbMatrix-----------------------------------------------------------------------')
    
    print('Input:')
    print('------')
    
    Q13 = ml.matrix("[0.1 0.5 0.4; 0.9 0.1 0; 0.3 -0.1 0.4]")
    print(Q13)    
    
    print('Test:')
    print('-----')
    
    print('CheckProbMatrix(Q13):')
    flag=CheckProbMatrix(Q13)
    print(flag)
    assert flag==0,'CheckProbMatrix did not recognize the negative entry!'
    
    print('Input:')
    print('------')
    
    Q14 = ml.matrix("[0.1 0.5 0.4; 0.9 0.1 0; 0.3 0.1 0.4]")
    print(Q14)    

    print('Test:')
    print('-----')
    
    print('CheckProbMatrix(Q14):')
    flag=CheckProbMatrix(Q14)
    print(flag)
    assert flag==0,'CheckProbMatrix did not recognize the invalid row sum!'
    
    print('Input:')
    print('------')
    
    Q15 = ml.matrix("[0.1 0.5 0.4; 0.9 0.1 0; 0.3 0.3 0.4]")
    print(Q15)    
    
    print('Test:')
    print('-----')
    
    print('CheckProbMatrix(Q15):')
    flag=CheckProbMatrix(Q15)
    print(flag)
    assert flag==1,'CheckProbMatrix did not recognize that the input is valid!'
    
    print('Input:')
    print('------')
    
    Q16 = ml.matrix("[0.1 0.5 0.4; 0.9 0.1 0; 0.3 0.3 0.4]")
    print(Q16)    
    
    print('Test:')
    print('-----')
    
    print('CheckProbMatrix(Q16,True):')
    flag=CheckProbMatrix(Q16,True)
    print(flag)
    assert flag==0,'CheckProbMatrix did not recognize wrong transient matrix!'
    
    print('Input:')
    print('------')
    
    Q17 = ml.matrix("[0.1 0.5 0.4; 0.9 0.1 0; 0.3 0.1 0.4]")
    print(Q17)    
    
    print('Test:')
    print('-----')
    
    print('CheckProbMatrix(Q17,True):')
    flag=CheckProbMatrix(Q17,True)
    print(flag)
    assert flag==1,'CheckProbMatrix did not recognize that the input is valid!'
    
    print('--CheckProbVector-----------------------------------------------------------------------')
    
    print('Input:')
    print('------')
    
    Q18 = ml.matrix("[1.1 -0.1]")
    print(Q18)
    
    print('Test:')
    print('-----')
    
    print('CheckProbVector(Q18):')
    flag=CheckProbVector(Q18)
    print(flag)
    assert flag==0,'CheckProbVector did not recognize the negative entry!'
    
    print('Input:')
    print('------')
    
    Q19 = ml.matrix("[1.1 0.1]")
    print(Q19)
    
    print('Test:')
    print('-----')
    
    print('CheckProbVector(Q19):')
    flag=CheckProbVector(Q19)
    print(flag)
    assert flag==0,'CheckProbVector did not recognize invalid sum!'
    
    print('Input:')
    print('------')
    
    Q20 = ml.matrix("[1 0]")
    print(Q20)
    
    print('Test:')
    print('-----')
    
    print('CheckProbVector(Q20):')
    flag=CheckProbVector(Q20)
    print(flag)
    assert flag==1,'CheckProbVector did not recognize that the input is valid!'
    
    print('Input:')
    print('------')
    
    Q21 = ml.matrix("[0.9 -0.1]")
    print(Q21)
    
    print('Test:')
    print('-----')
    
    print('CheckProbVector(Q21,True):')
    flag=CheckProbVector(Q21,True)
    print(flag)
    assert flag==0,'CheckProbVector did not recognize the negative entry!'
    
    print('Input:')
    print('------')
    
    Q22 = ml.matrix("[0.9 0.1]")
    print(Q22)
    
    print('Test:')
    print('-----')
    
    print('CheckProbVector(Q22,True):')
    flag=CheckProbVector(Q22,True)
    print(flag)
    assert flag==1,'CheckProbVector did not recognize that the prob. vector is not transient!'
    
    print('Input:')
    print('------')
    
    Q23 = ml.matrix("[0.8 0.1]")
    print(Q23)
    
    print('Test:')
    print('-----')
    
    print('CheckProbVector(Q23,True):')
    flag=CheckProbVector(Q23,True)
    print(flag)
    assert flag==1,'CheckProbVector did not recognize that the input is valid!'
