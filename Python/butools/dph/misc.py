# -*- coding: utf-8 -*-
"""
Created on Sun Aug 24 15:47:28 2014

@author: gabor
"""

import subprocess
import os
import os.path
import butools
from butools.dph import CheckDPHRepresentation
import numpy as np
import numpy.matlib as ml
from numpy.random import rand

def SamplesFromDPH (a,A,k):
    """
    Generates random samples from a discrete phase-type 
    distribution.
    
    Parameters
    ----------
    alpha : matrix, shape (1,M)
        The initial probability vector of the discrete phase-
        type distribution.
    A : matrix, shape (M,M)
        The transition probability  matrix of the discrete phase-
        type distribution.
    K : integer
        The number of samples to generate.
    prec : double, optional
        Numerical precision to check if the input phase-type
        distribution is valid. The default value is 1e-14.
    
    Returns
    -------
    x : vector, length(K)
        The vector of random samples
    """

    if butools.checkInput and not CheckDPHRepresentation(a,A):
        raise Exception("SamplesFromDPH: input is not a valid DPH representation!")

    # auxilary variables
    a = a.A.flatten()
    N = len(a)
    cummInitial = np.cumsum(a)
    logp = np.log(np.diag(A));
    sojourn = 1.0/(1.0-np.diag(A))
    nextpr = ml.matrix(np.diag(sojourn))*A
    nextpr = nextpr - ml.matrix(np.diag(np.diag(nextpr)))
    nextpr = np.hstack((nextpr, 1.0-np.sum(nextpr,1)))
    nextpr = np.cumsum(nextpr,1)
    
    x = np.zeros(k)
    for n in range(k):
        time = 0

        # draw initial distribution
        r = rand()
        state = 0
        while cummInitial[state]<=r:
            state+=1

        # play state transitions
        while state<N:
            time += 1 + int(np.log(rand()) / logp[state])
            r = rand()
            nstate = 0
            while nextpr[state,nstate]<=r:
                nstate += 1
            state = nstate
        x[n] = time
    return x
    
def ImageFromDPH(alpha,A,outFileName=None,prec=1e-13):
    """
    Depicts the given discrete phase-type distribution,
    and either displays it or saves it to file.
    
    Parameters
    ----------
    alpha : matrix, shape (1,M)
        The initial probability vector of the discrete phase-
        type distribution.
    A : matrix, shape (M,M)
        The transition probability  matrix of the discrete phase-
        type distribution.
    outFileName : string, optional
        If it is not provided, or equals to 'display', the
        image is displayed on the screen, otherwise it is 
        written to the file. The file format is deduced 
        from the file name.
    prec : double, optional
        Transition probabilities less then prec are 
        considered to be zero and are left out from the 
        image. The default value is 1e-13.
    
    Notes
    -----
    The 'graphviz' software must be installed and available
    in the path to use this feature.
    """

    if butools.checkInput and not CheckDPHRepresentation(alpha,A):
        raise Exception("ImageFromDPH: input is not a valid DPH representation!")

    if outFileName==None or outFileName=="display":
        outputFile = ".result.png"
        displ = True
    else:
        outputFile = outFileName
        displ = False
    
    inputFile = "_temp.dot"

    fid = open(inputFile,'w')
    fid.write ("digraph G {\n")
    fid.write ("\trankdir=LR;\n")
    fid.write ('\tnode [shape=circle,width=0.3,height=0.3,label=""];\n')
    
    # nodes
    alpha = alpha.A.flatten()
    for i in range(len(alpha)):
        fid.write ('\tn{0} [xlabel=<<i>{1}</i>>];\n'.format(i, alpha[i]))
    
    # transitions to a non-absorbing state
    for i in range(len(alpha)):
        for j in range(len(alpha)):
            if abs(A[i,j])>prec:
                fid.write('\tn{0} -> n{1} [label="{2}"];\n'.format(i, j, A[i,j]))
    
    # transitions to the absorbing state
    fid.write('\tab [style=filled];\n')
    a = 1.0-np.sum(A,1).A.flatten()
    for i in range(len(alpha)):
        if abs(a[i])>prec:
            fid.write('\tn{0} -> ab [label="{1}"];\n'.format(i, a[i]))

    fid.write('}\n')
    fid.close()

    ext = os.path.splitext(outputFile)[1]
    subprocess.call(['dot', "-T"+ext[1:], inputFile, "-o", outputFile])
    os.remove (inputFile)
   
    if displ:
        from IPython.display import Image
        i = Image(filename=outputFile)
        os.remove(outputFile)        
        return i
