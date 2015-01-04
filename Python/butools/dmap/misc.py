# -*- coding: utf-8 -*-
"""
Created on Sun Aug 31 17:33:24 2014

@author: gabor
"""

import butools
import numpy as np
import numpy.matlib as ml
from butools.dmap import CheckDMAPRepresentation, CheckDMMAPRepresentation
from numpy.random import rand
from butools.mc import DTMCSolve
from butools.utils import SumMatrixList


def SamplesFromDMMAP (D, k, initial=None, prec=1e-14):

    if butools.checkInput and not CheckDMMAPRepresentation (D, prec):
        raise Exception("SamplesFromDMMAP: Input is not a valid DMMAP representation!")    

    N = D[0].shape[0]
    
    if initial==None:
        # draw initial state according to the stationary distribution
        stst = DTMCSolve(SumMatrixList(D)).A.flatten()
        cummInitial = np.cumsum(stst)
        r = rand()
        state = 0
        while cummInitial[state]<=r:
            state+=1
    else:
        state = initial

    # auxilary variables
    sojourn = 1.0/(1.0-np.diag(D[0]))
    logp = np.log(np.diag(D[0]))
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
            time += 1 + int(np.log(rand()) / logp[state])
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

def SamplesFromDMAP (D0, D1, k, initial=None, prec=1e-14):

    if butools.checkInput and not CheckDMAPRepresentation (D0, D1, prec):
        raise Exception("SamplesFromDMAP: Input is not a valid DMAP representation!")    

    return SamplesFromDMMAP((D0,D1),k,initial,prec);

import os
import os.path
from subprocess import call

def ImageFromDMMAP (D, outFileName="display", prec=1e-13):

    if outFileName=="display":
        outputFile = ".result.png"
        displ = True
    else:
        outputFile = outFileName
        displ = False
    
    inputFile = "_temp.dot"

    f = open(inputFile,"w")
    f.write("digraph G {\n")
    f.write("\trankdir=LR;\n")
    f.write('\tnode [shape=circle,width=0.3,height=0.3,label=""];\n')

    N = D[0].shape[0]
    
    # transitions without arrivals
    Dx=D[0]
    for i in range(N):
        for j in range(N):
            if abs(Dx[i,j])>prec:
                f.write('\tn{0} -> n{1} [label="{2}"];\n'.format(i, j, Dx[i,j]))
    
    # transitions with arrivals
    for k in range(1,len(D)):
        Dx=D[k]
        for i in range(N):
            for j in range(N):
                if abs(Dx[i,j])>prec:
                    if len(D)==2:
                        f.write('\tn{0} -> n{1} [style="dashed",label="{2}"];\n'.format(i, j, Dx[i,j]))
                    else:
                        f.write('\tn{0} -> n{1} [style="solid",fontcolor="/dark28/{2}",color="/dark28/{3}",label="{4}"];\n'.format(i, j, min(k-1,8), min(k-1,8), Dx[i,j]))

    f.write('}\n')
    f.close()

    ext = os.path.splitext(outputFile)[1]
    call(['dot', '-T' + ext[1:], "_temp.dot", '-o', outputFile])   
    os.remove (inputFile)
   
    if displ:
        from IPython.display import Image
        i = Image(filename=outputFile)
        os.remove(outputFile)        
        return i

def ImageFromDMAP (D0, D1, outFileName="display", prec=1e-13):

    return ImageFromDMMAP((D0,D1),outFileName,prec)