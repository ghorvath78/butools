# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 17:33:10 2014

@author: gabor
"""

import math
import scipy.linalg as la
import butools

def CheckMoments (m, prec=1e-14):
    """
    Checks if the given moment sequence is valid in the sense
    that it belongs to a distribution with support (0,inf).
    
    This procedure checks the determinant of `\Delta_n`
    and `\Delta_n^{(1)}` according to [1]_.
    
    Parameters
    ----------
    m : list of doubles, length 2N+1
        The (raw) moments to check 
        (starts with the first moment).
        Its length must be odd.
    prec : double, optional
        Entries with absolute value less than prec are 
        considered to be zeros. The default value is 1e-14.
        
    Returns
    -------
    r : bool
        The result of the check.
    
    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Stieltjes_moment_problem
    """

    if butools.checkInput and len(m)%2==0:
        raise Exception("CheckMoments: the number of moments must be odd!")
    
    m = [1.0] + m
    N = math.floor(len(m)/2)-1
    
    for n in range(N+1):
        H = la.hankel(m[0:n+1], m[n:2*n+1])
        H0 = la.hankel(m[1:n+2], m[n+1:2*n+2])
        if la.det(H)<-prec or la.det(H0)<-butools.checkPrecision:
            return False
    return True