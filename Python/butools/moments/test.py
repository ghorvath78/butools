# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 17:20:24 2014

@author: gabor
"""

import butools
import numpy as np
from butools.moments import *

def TestMomentsPackage ():


    print('---BuTools: Moments package test file---')

    print('Enable the verbose messages with the BuToolsVerbose flag')
    butools.verbose = True
    
    print('Enable input parameter checking with the BuToolsCheckInput flag')
    butools.checkInput = True


    print('--NormMomsFromMoms/MomsFromNormMoms-------------------------------------------------')

    print('Test:')
    print('-----')
    M = [1.2, 5, 38, 495, 9215]
    print(M)
    print('nmoms=NormMomsFromMoms(M)')
    nmoms=NormMomsFromMoms(M)
    print(nmoms)
    print('moms=MomsFromNormMoms(nmoms)')
    moms=MomsFromNormMoms(nmoms)
    print(moms)
    assert np.max(np.abs(np.array(moms)-np.array(M)))<1e-11, 'Calling the moment conversion and its inverse did not give back the original moments!'

    print('--ReducedMomsFromMoms/MomsFromReducedMoms-------------------------------------------')

    print('Test:')
    print('-----')
    print(M)
    print('rmoms=ReducedMomsFromMoms(M)')
    rmoms=ReducedMomsFromMoms(M)
    print(rmoms)
    print('moms=MomsFromReducedMoms(rmoms)')
    moms=MomsFromReducedMoms(rmoms)
    print(moms)
    assert np.max(np.abs(np.array(moms)-np.array(M)))<1e-11, 'Calling the moment conversion and its inverse did not give back the original moments!'

    print('--FactorialMomsFromMoms/MomsFromFactorialMoms---------------------------------------')

    print('Test:')
    print('-----')
    M = [1.3, 2.4, 6.03, 20.5, 89.5, 474.9]
    print(M)
    print('fmoms=FactorialMomsFromMoms(M)')
    fmoms=FactorialMomsFromMoms(M)
    print(fmoms)
    print('moms=MomsFromFactorialmoms(fmoms)')
    moms=MomsFromFactorialMoms(fmoms)
    print(moms)
    assert np.max(np.abs(np.array(moms)-np.array(M)))<1e-11, 'Calling the moment conversion and its inverse did not give back the original moments!'


    print('--JFactorialMomsFromJMoms/JMomsFromJFactorialMoms-----------------------------------')

    print('Test:')
    print('-----')
    N = [[0.7, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]]
    print(N)
    print('JFmoms=JFactorialMomsFromJMoms(N)')
    JFmoms=JFactorialMomsFromJMoms(N)
    print(JFmoms)
    print('Jmoms=JMomsFromJFactorialMoms(JFmoms)')
    Jmoms=JMomsFromJFactorialMoms(JFmoms)
    print(Jmoms)
    assert np.max(np.max(np.abs(np.array(moms)-np.array(M))))<1e-11, 'Calling the moment conversion and its inverse did not give back the original moments!'

    print('--HankelMomsFromMoms/MomsFromHankelMoms---------------------------------------------')

    print('Test:')
    print('-----')
    print(M)
    print('hmoms=HankelMomsFromMoms(M)')
    hmoms=HankelMomsFromMoms(M)
    print(hmoms)
    print('moms=MomsFromHankelMoms(hmoms)')
    moms=MomsFromHankelMoms(hmoms)
    print(moms)
    assert np.max(np.abs(np.array(moms)-np.array(M)))<1e-11, 'Calling the moment conversion and its inverse did not give back the original moments!'

    print('--CheckMoments----------------------------------------------------------------------')

    print('Test:')
    print('-----')
    M = [1.2, 5, 8, 29, 3412]
    print(M)
    print('flag=CheckMoments(M)')
    flag=CheckMoments(M)
    print(flag)
    assert flag==False, 'CheckMoments did not recognize a valid moment sequence!'
    M = [1.3, 2.4, 6.03, 20.5, 89.5]
    print(M)
    print('flag=CheckMoments(M)')
    flag=CheckMoments(M)
    print(flag)
    assert flag==True, 'CheckMoments did not recognize an invalid moment sequence!'

if __name__ == "__main__":
    TestMomentsPackage()