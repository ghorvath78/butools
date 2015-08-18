# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 08:24:03 2013

@author: Gabor Horvath
"""
import numpy as np
import numpy.matlib as ml
import math
import butools
import sys
from IPython.display import clear_output

def PHFromTrace (trace, orders, weights=[], maxIter=200, stopCond=1e-7, initial=None, result="vecmat", retlogli=True):
    """
    Performs PH distribution fitting using the EM algorithm
    (G-FIT, [1]_).
    
    Parameters
    ----------
    trace : column vector, length K
        The samples of the trace
    orders : list of int, length(N), or int
        The length of the list determines the number of 
        Erlang branches to use in the fitting method.
        The entries of the list are the orders of the 
        Erlang distributions. If this parameter is a 
        single integer, all possible branch number - order
        combinations are tested where the total number of 
        states is "orders".
    maxIter : int, optional
        Maximum number of iterations. The default value is
        200
    stopCond : double, optional
        The algorithm stops if the relative improvement of
        the log likelihood falls below stopCond. The 
        default value is 1e-7
    initial : tuple of two vectors, optional
        The initial values of the branch probabilities and
        rate parameters is given by this tuple. If not 
        given, a default initial guess is determined and 
        the algorithm starts from there.
    result : {"vecmat", "vecvec"}, optional
        The result can be returned two ways. If "vecmat" is
        selected, the result is returned in the classical
        representation of phase-type distributions, thus the
        initial vector and the generator matrix. 
        If "vecvec" is selected, two vectors are returned, 
        one holds the branch probabilities, and the second
        holds the rate parameters of the Erlang branches.
        The default value is "vecmat"
    
    Returns
    -------
    (alpha, A) : tuple of matrix, shape (1,M) and matrix, shape (M,M)
        If the "vecmat" result format is chosen, the function
        returns the initial probability vector and the
        generator matrix of the phase type distribution.
    (pi, lambda) : tuple of vector, length N and vector, length N
        If the "vecvec" result format is chosen, the function
        returns the vector of branch probabilities and the
        vector of branch rates in a tuple.
    logli : double
        The log-likelihood divided by the trace length
        
    Notes
    -----
    This procedure is quite fast in the supported 
    mathematical frameworks. If the maximum speed is
    needed, please use the multi-core optimized c++
    implementation called SPEM-FIT_.
    
    .. _SPEM-FIT: https://bitbucket.org/ghorvath78/spemfit
    
    References
    ----------
    .. [1] Thummler, Axel, Peter Buchholz, and Miklós Telek.
           A novel approach for fitting probability 
           distributions to real trace data with the EM 
           algorithm. Dependable Systems and Networks, 2005.
    """
    
    def allorders (branches, sumorders):
        if branches==1:
            return [[sumorders]]
        else:
            o = [];
            for i in range(sumorders-branches+1):
                x = allorders (branches-1, sumorders-i-1)
                for j in range(len(x)):
                    xt = x[j]
                    xt.append(i+1)
                    xt.sort()
                    # check if we have it already
                    if o.count(xt)==0:
                        o.append(xt)
#                    for ok in o:
#                        if ok==xt
#                            break;
#                    else:
#                        o.append(xt)
            return o
               
    
    if type(orders) is int:
        bestres = ([],[],-np.inf)
        for br in range(2,orders+1):
            allord = allorders(br, orders)
            for ordk in allord:
                res = PHFromTrace (trace, ordk, weights, maxIter, stopCond, initial, result, True)
                if res[2] > bestres[2]:
                    bestres = res
        if retlogli:
            return bestres
        else:
            return (bestres[0], bestres[1])
    
    M = len(orders)
    K = len(trace)

    # initial alpha and lambda is such that the mean is matched
    if initial==None:
        alphav = np.ones(M) / M
        lambd = orders * np.linspace(1,M,M)
        trm = np.sum(trace)/len(trace)
        inim = np.sum(alphav / np.linspace(1,M,M))
        lambd = lambd * inim / trm
    elif len(initial)==2:
        if len(initial[0])==M and len(initial[1])==M:
            alphav = initial[0]
            lambd = initial[1]
        else:
            raise Exception("The length of the initial branch probability and rate vectors is not consistent with the length of the orders vector!")
    else:
        raise Exception("Invalid initial branch probability and rate vectors!")

    if weights==[]:
        weights = np.ones(len(trace))/K
    W = weights / np.sum(weights)

    Q = np.zeros((M, K))
    logli, ologli = 1e-14, 1
    steps = 1
    while abs((ologli-logli)/logli) > stopCond and steps<=maxIter:
        ologli = logli
        # E-step:
        for i in range (M):
            Q[i,:] = (alphav[i]*(lambd[i]*trace)**(orders[i]-1) / math.factorial(orders[i]-1) * lambd[i]) * np.exp(-lambd[i]*trace)
        nor = np.sum(Q,0)
        for i in range (M):
            Q[i,:] /= nor
        logli = np.log(nor).dot(W) #np.sum(np.log(nor)) / K
        if butools.verbose and steps%10==0:
            clear_output ()
            print ("iteration: ", steps, ", logli: ", logli)
            sys.stdout.flush()
        # M-step:
        v1 = Q.dot(W) # np.sum(Q,1)
        v2 = Q.dot(trace*W) # Q.dot(trace)
        alphav = v1 # v1/K
        lambd = orders*v1 / v2
        steps += 1

    if butools.verbose:
        clear_output ()
        print ("EM algorithm terminated.", orders)
        print ("Num of iterations: ", steps, ", logli: ", logli)
        sys.stdout.flush()

    if result=="vecvec":
        if retlogli:
            return (alphav, lambd, logli)
        else:
            return (alphav, lambd)
    elif result=="vecmat":
        # construct the vector and the matrix representation
        N = sum(orders)
        alpha = ml.zeros((1,N))
        A = ml.zeros ((N,N))
        ix = 0
        for i in range(M):
            alpha[0,ix] = alphav[i]
            if orders[i]==1:
                A[ix,ix] = -lambd[i]
            else:
                A[ix:ix+orders[i], ix:ix+orders[i]] = lambd[i]*(np.diag(np.ones(orders[i]-1),1)-np.diag(np.ones(orders[i])))
            ix += orders[i]       
        if retlogli:
            return (alpha, A, logli)
        else:
            return (alpha, A)
    else:
        raise Exception("Unknown result format given! (valid are: vecmat and vecvec)")
