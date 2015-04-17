# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 08:24:03 2013

@author: Gabor Horvath
"""
import numpy as np
import math

def SquaredDifference (p1, p2):
    """
    Returns the squared difference between two vectors.
    
    Parameters
    ----------
    p1 : vector, length M
        The first vector
    p2 : vector, length M
        The second vector
    
    Returns
    -------
    sd : double
        The squared difference calculated as
        `sq=\sum_{i=1}^M (p1_i-p2_i)^2`
    """

    return np.sum(np.square(np.array(p1)-np.array(p2)))

def EmpiricalSquaredDifference (f1, f2, intBounds):
    """
    Returns the squared difference of two continuous 
    functions given by samples and the bounds of the 
    corresponding intervalls.
    
    This function can be used to characterize the distance
    between two density functions, distribution functions, 
    etc.
    
    Parameters
    ----------
    f1 : vector, length M
        Samples of the first continuous function
    f2 : vector, length M
        Samples of the second continuous function
    intBounds : vector, length M+1
        The bounds of the intervals. The ith sample
        corresponds to the interval 
        (intbounds(i),intbounds(i+1))
    
    Returns
    -------
    sd : double
        The squared difference
    """

    intlens = intBounds[1:] - intBounds[0:-1]
    p1 = f1 * intlens
    p2 = f2 * intlens
    return SquaredDifference (p1, p2)
    
def RelativeEntropy (p1, p2):
    """
    Returns the relative entropy (aka Kullback–Leibler
    divergence) of two vectors.
    
    Parameters
    ----------
    p1 : vector, length M
        The first vector
    p2 : vector, length M
        The second vector
    
    Returns
    -------
    re : double
        The relative entropy calculated as
        `re=\sum_{i=1}^M p1_i |\log(p1_i/p2_i)|`
    """

    re = 0
    for i in range(len(p1)):
        if p1[i] > 0.0:
            re += p1[i]*abs(math.log(p1[i]/p2[i]))
    return re
#    return np.sum(p1*np.log(p1/p2))
    
def EmpiricalRelativeEntropy (f1, f2, intBounds):
    """
    Returns the relative entropy (aka Kullback–Leibler
    divergence) of two continuous functions given by samples
    and the bounds of the corresponding intervalls.
    
    This function can be used to characterize the distance
    between two density functions, distribution functions, 
    etc.
    
    Parameters
    ----------
    f1 : vector, length M
        Samples of the first continuous function
    f2 : vector, length M
        Samples of the second continuous function
    intBounds : vector, length M+1
        The bounds of the intervals. The ith sample
        corresponds to the interval 
        (intbounds(i),intbounds(i+1))
    
    Returns
    -------
    re : double
        The relative entropy
    """

    intlens = intBounds[1:] - intBounds[0:-1]
    p1 = f1 * intlens
    p2 = f2 * intlens
    return RelativeEntropy (p1, p2)

