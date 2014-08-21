# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 17:33:10 2014

@author: gabor
"""

import math
import scipy.linalg as la
import butools

def CheckMoments (m, prec=1e-14):

    if butools.checkInput and len(m)%2==0:
        raise Exception("CheckMoments: the number of moments must be odd!")
    
    m = [1.0] + m
    N = math.floor(len(m)/2)-1
    
    for n in range(N+1):
        H = la.hankel(m[0:n+1], m[n:2*n+1])
        H0 = la.hankel(m[1:n+2], m[n+1:2*n+2])
        if la.det(H)<-prec or la.det(H0)<-prec:
            return False
    return True