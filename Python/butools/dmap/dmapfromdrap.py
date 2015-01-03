# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 08:21:05 2013

@author: gabor
"""

import numpy as np
from numpy import linalg as la
import numpy.matlib as ml
import butools
from butools.reptrans import FindMarkovianRepresentation
from butools.dmap import CheckDRAPRepresentation, CheckDMRAPRepresentation

def DMMAPFromDMRAP (H, prec=1e-14):

    if butools.checkInput and not CheckDMRAPRepresentation (H, prec):
        raise Exception("DMMAPFromDMRAP: Input is not a valid DMRAP representation!")    

    def transfun (oH, B):
        return [la.inv(B)*oHk*B for oHk in oH]
        
    def evalfun (oH, k=0):
        Ones = ml.ones(oH[0].shape)
        if k%2 == 0:
            dist = np.min(oH[0])
            for oHk in oH:
                dist = min(dist, np.min(oHk), np.min(Ones-oHk))
            return -dist
        else:
            dist = np.sum(oH[0][oH[0]<0])
            for oHk in oH:
                dist += min(np.sum(oHk[oHk<0]), np.sum(oHk[Ones-oHk<0]))
            return -dist

    return FindMarkovianRepresentation (H, transfun, evalfun, prec)

def DMAPFromDRAP (H0, H1, prec=1e-14):

    if butools.checkInput and not CheckDRAPRepresentation (H0, H1, prec):
        raise Exception("DMAPFromDRAP: Input is not a valid DRAP representation!")    

    return DMMAPFromDMRAP([H0,H1], prec)
