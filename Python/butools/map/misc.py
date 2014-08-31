# -*- coding: utf-8 -*-
"""
Created on Sun Aug 31 17:33:24 2014

@author: gabor
"""

import butools
import numpy as np
import numpy.matlib as ml
from butools.map import CheckMAPRepresentation, CheckMMAPRepresentation, MarginalDistributionFromMMAP
from numpy.random import rand


def SamplesFromMMAP (D, k, initial=None, prec=1e-14):

    if butools.checkInput and not CheckMMAPRepresentation (D, prec):
        raise Exception("SamplesFromMMAP: Input is not a valid MMAP representation!")    

    N = D[0].shape[0]
    
    if initial==None:
        # draw initial state according to the stationary distribution
        stst = MarginalDistributionFromMMAP(D).A.flatten()
        cummInitial = np.cumsum(stst)
        r = rand()
        state = 0
        while cummInitial[state]<=r:
            state+=1
    else:
        state = initial

    # auxilary variables
    sojourn = -1.0/np.diag(D[0])
    nextpr = ml.matrix(np.diag(sojourn))*D[0]
    nextpr = nextpr - ml.matrix(np.diag(np.diag(nextpr)))
    for i in range(1,len(D)):
        nextpr = np.hstack((nextpr, np.diag(sojourn)*D[i]))
    nextpr = np.cumsum(nextpr,1)
    
    if len(D)>2:
        x = np.empty((k,2))
    else:
        x = np.empty(k)

    for n in range(k):
        time = 0

        # play state transitions
        while state<N :
            time -= np.log(rand()) * sojourn[state]
            r = rand()
            nstate = 0
            while nextpr[state,nstate]<=r:
                nstate += 1
            state = nstate
        if len(D)>2:
            x[n,0] = time
            x[n,1] = state//N
        else:
            x[n] = time
        state = state % N
    
    return x

def SamplesFromMAP (D0, D1, k, initial=None, prec=1e-14):

    if butools.checkInput and not CheckMAPRepresentation (D0, D1, prec):
        raise Exception("SamplesFromMAP: Input is not a valid MAP representation!")    

    return SamplesFromMMAP((D0,D1),k,initial,prec);
