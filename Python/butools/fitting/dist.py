# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 08:24:03 2013

@author: Gabor Horvath
"""
import numpy as np
import math

def SquaredDifference (p1, p2):

    return np.sum(np.square(np.array(p1)-np.array(p2)))

def EmpiricalSquaredDifference (f1, f2, intBounds):

    intlens = intBounds[1:] - intBounds[0:-1]
    p1 = f1 * intlens
    p2 = f2 * intlens
    return SquaredDifference (p1, p2)
    
def RelativeEntropy (p1, p2):

    re = 0
    for i in range(len(p1)):
        if p1[i] > 0.0:
            re += p1[i]*abs(math.log(p1[i]/p2[i]))
    return re
#    return np.sum(p1*np.log(p1/p2))
    
def EmpiricalRelativeEntropy (f1, f2, intBounds):

    intlens = intBounds[1:] - intBounds[0:-1]
    p1 = f1 * intlens
    p2 = f2 * intlens
    return RelativeEntropy (p1, p2)

