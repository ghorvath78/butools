# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 08:24:03 2013

@author: Gabor Horvath
"""
import numpy as np
import numpy.matlib as ml
import math
import butools
from butools.mc import DTMCSolve
import sys
from IPython.display import clear_output

def MAPFromTrace (trace, orders, maxIter=200, stopCond=1e-7, initial=None, result="matmat", retlogli=True):
    """
    Performs MAP fitting using the EM algorithm (ErCHMM, 
    [1]_, [2]_).
    
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
    initial : tuple of a vector and a matrix, shape(N,N), optional
        The rate parameters of the Erlang distributions 
        and the branch transition probability matrix to be
        used initially. If not given, a default initial 
        guess is determined and the algorithm starts from 
        there.
    result : {"vecmat", "matmat"}, optional
        The result can be returned two ways. If "matmat" is
        selected, the result is returned in the classical
        representation of MAPs, thus the D0 and D1 matrices.
        If "vecmat" is selected, the rate parameters of the
        Erlang branches and the branch transition probability
        matrix are returned. The default value is "matmat"
    
    Returns
    -------
    (D0, D1) : tuple of matrix, shape (M,M) and matrix, shape (M,M)
        If the "matmat" result format is chosen, the function
        returns the D0 and D1 matrices of the MAP
    (lambda, P) : tuple of vector, length N and matrix, shape (M,M)
        If the "vecmat" result format is chosen, the function
        returns the vector of the Erlang rate parameters of 
        the branches and the branch transition probability 
        matrix
    logli : double
        The log-likelihood divided by the trace length
        
    Notes
    -----
    This procedure is quite slow in the supported 
    mathematical frameworks. If the maximum speed is
    needed, please use the multi-core optimized c++
    implementation called SPEM-FIT_.
    
    .. _SPEM-FIT: https://bitbucket.org/ghorvath78/spemfit
    
    References
    ----------
    .. [1] Okamura, Hiroyuki, and Tadashi Dohi. Faster 
           maximum likelihood estimation algorithms for 
           Markovian arrival processes. Quantitative 
           Evaluation of Systems, 2009. QEST'09. Sixth 
           International Conference on the. IEEE, 2009.
    
    .. [2] Horváth, Gábor, and Hiroyuki Okamura. A Fast EM
           Algorithm for Fitting Marked Markovian Arrival 
           Processes with a New Special Structure. Computer
           Performance Engineering. Springer Berlin 
           Heidelberg, 2013. 119-133.
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
        bestOrders = []
        for br in range(2,orders+1):
            allord = allorders(br, orders)
            for ordk in allord:
                if butools.verbose:
                    print("Trying orders ", ordk)
                res = MAPFromTrace (trace, ordk, maxIter, stopCond, initial, result, True)
                if res[2] > bestres[2]:
                    bestres = res
                    bestOrders = ordk
        if butools.verbose:
            print("Best solution: logli =", bestres[2],"orders =", bestOrders)
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
        P = (ml.ones((M,1))*ml.matrix([alphav])).A
    elif len(initial)==2:
        if len(initial[0])==M and initial[1].shape==(M,M):
            lambd = initial[0]
            P = np.array(initial[1])
            alphav = DTMCSolve(ml.matrix(P))
        else:
            raise Exception("The length of the initial branch probability and rate vectors is not consistent with the length of the orders vector!")
    else:
        raise Exception("Invalid initial branch probability and rate vectors!")

    Q = np.zeros((M, K))
    A = np.zeros((K, M))
    B = np.zeros((M, K))
    Ascale = np.zeros(K)
    Bscale = np.zeros(K)
    logli, ologli = 1e-14, 0
    steps = 1
    while abs((ologli-logli)/logli) > stopCond and steps<=maxIter:
        ologli = logli
        # E-step:
        for i in range (M):
            Q[i,:] = ((lambd[i]*trace)**(orders[i]-1) / math.factorial(orders[i]-1) * lambd[i]) * np.exp(-lambd[i]*trace)
        # forward likelihood vectors:
        prev = alphav
        scprev = 0
        for k in range(K):
            prev = prev.dot(np.diag(Q[:,k])).dot(P)
            scale = math.log2(np.sum(prev))
            prev = prev * 2**-scale
            Ascale[k] = scprev + scale
            A[k,:] = prev
            scprev = Ascale[k]
        Av = np.vstack((alphav, A[0:-1,:]))
        Ascalev = np.hstack(([0], Ascale[0:-1]))
        # backward likelihood vectors:
        nnext = np.ones(M)
        scprev = 0
        for k in range(K-1,-1,-1):
            nnext = np.diag(Q[:,k]).dot(P).dot(nnext)
            scale = math.log2(np.sum(nnext))
            nnext = nnext * 2**-scale
            Bscale[k] = scprev + scale
            B[:,k] = nnext
            scprev = Bscale[k]
        Bv = np.hstack((B[:,1:], np.ones((M,1))))
        Bscalev = np.hstack((Bscale[1:], [0]))

        llh = alphav.dot(B[:,0])
        logli = (math.log(llh) + Bscale[0] * math.log(2)) / K
        illh = 1.0 / llh

        # M-step:
        # Calculate new estimates for the parameters
        AB = Av*B.T;
        nor = np.sum(AB,1)
        for m in range(M):
            AB[:,m] /= nor
        v1 = np.sum(AB,0)
        v2 = AB.T.dot(trace).T
        alphav = v1/K
        lambd = (orders*v1 / v2).T
        
        Avv = Av*Q.T
        nor = illh*2**(Ascalev+Bscalev-Bscale[0]).T
        for m in range(M):
            Avv[:,m] *= nor
        P = (Avv.T.dot(Bv.T))*P
        for m in range(M):
            P[m,:] /= np.sum(P[m,:])

        steps += 1
        if butools.verbose and steps%10==0:
            clear_output ()
            print ("iteration: ", steps, ", logli: ", logli)
            sys.stdout.flush()

    if butools.verbose:
        clear_output ()
        print ("Num of iterations: ", steps, ", logli: ", logli)
        print ("EM algorithm terminated.", orders)
        sys.stdout.flush()

    if result=="vecmat":
        if retlogli:
            return (lambd, P, logli)
        else:
            return (lambd, P)
    elif result=="matmat":
        N = np.sum(orders)
        D0 = ml.zeros ((N,N))
        ix = 0
        for i in range(M):
            if orders[i]==1:
                D0[ix,ix] = -lambd[i]
            else:
                D0[ix:ix+orders[i], ix:ix+orders[i]] = lambd[i]*(np.diag(np.ones(orders[i]-1),1)-np.diag(np.ones(orders[i])))
            ix += orders[i]       

        D1 = ml.zeros((N,N))
        indicesTo = np.hstack(([0], np.cumsum(orders[0:-1])))
        indicesFrom = np.cumsum(orders)-1
        X = np.diag(lambd).dot(P)
        ix = 0
        for i in indicesFrom:
            jx = 0
            for j in indicesTo:
                D1[i,j] = X[ix, jx]
                jx += 1
            ix += 1
        if retlogli:
            return (D0, D1, logli)
        else:
            return (D0, D1)
    else:
        raise Exception("Unknown result format given! (valid are: vecmat and vecvec)")
