# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 11:01:03 2013

@author: gabor
"""

import numpy.matlib as ml
import numpy.linalg as la
import butools
from butools.map.basemap import MarginalDistributionFromMRAP
from butools.reptrans import MStaircase
from butools.map import CheckRAPRepresentation, CheckMRAPRepresentation

def MinimalRepFromMRAP (H, how="obscont", precision=1e-12):

    if butools.checkInput and not CheckMRAPRepresentation (H, precision):
        raise Exception("MinimalRepFromMRAP: Input is not a valid MRAP representation!")    

    if how=="cont":
        B, n = MStaircase (H, ml.ones((H[0].shape[0],1)), precision)
        return [(la.inv(B)*Hk*B)[0:n,0:n] for Hk in H]
    elif how=="obs":
        alpha, A = MarginalDistributionFromMRAP (H)
        G = [Hk.T for Hk in H]
        B, n = MStaircase (G, alpha.T, precision)
        return [(la.inv(B)*Hk*B)[0:n,0:n] for Hk in H]
    elif how=="obscont":
        D = MinimalRepFromMRAP(H, "cont", precision)
        return MinimalRepFromMRAP(D, "obs", precision)
   
def MinimalRepFromRAP (H0, H1, how="obscont", precision=1e-12):

    if butools.checkInput and not CheckRAPRepresentation (H0, H1, precision):
        raise Exception("MinimalRepFromRAP: Input is not a valid RAP representation!")    

    return MinimalRepFromMRAP ([H0,H1], how, precision)
