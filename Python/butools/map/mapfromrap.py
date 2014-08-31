# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 08:21:05 2013

@author: gabor
"""

import numpy as np
from numpy import linalg as la
import butools
from butools.reptrans import FindMarkovianRepresentation
from butools.map import CheckRAPRepresentation, CheckMRAPRepresentation

def MMAPFromMRAP (H, prec=1e-14):

    if butools.checkInput and not CheckMRAPRepresentation (H, prec):
        raise Exception("MMAPFromMRAP: Input is not a valid MRAP representation!")    

    def transfun (oH, B):
        return [la.inv(B)*oHk*B for oHk in oH]
        
    def evalfun (oH, k=0):
        oH0 = oH[0] - np.diag(np.diag(oH[0]))
        if k%2 == 0:
            dist = np.min(oH0)
            for oHk in oH[1:]:
                dist = min(dist, np.min(oHk))
            return -dist
        else:
            dist = np.sum(oH0[oH0<0])
            for oHk in oH[1:]:
                dist += np.sum(oHk[oHk<0])
            return -dist

    return FindMarkovianRepresentation (H, transfun, evalfun, prec)

def MAPFromRAP (H0, H1, prec=1e-14):

    if butools.checkInput and not CheckRAPRepresentation (H0, H1, prec):
        raise Exception("MAPFromRAP: Input is not a valid RAP representation!")    

    return MMAPFromMRAP([H0,H1], prec)
