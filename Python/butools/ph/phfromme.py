# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 08:21:05 2013

@author: gabor
"""

import numpy as np
from numpy import linalg as la
from butools.reptrans import FindMarkovianRepresentation

def PHFromME (alpha, A, precision=1e-14):

    def transfun (orep, B):
        ao, Ao = orep
        return (ao*B, la.inv(B)*Ao*B)
        
    def evalfun (orep, k=0):
        ao, Ao = orep
        av= np.sum(-Ao,1)
        Ad = Ao-np.diag(np.diag(Ao))
        if k%2 == 0:
            return -min(np.min(ao), np.min(av), np.min(Ad))
        else:
            return -np.sum(ao[ao<0]) - np.sum(av[av<0]) - np.sum(Ad[Ad<0])

    return FindMarkovianRepresentation ((alpha,A), transfun, evalfun, precision)
