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

def PHFromTrace (trace, orders, weights=[], maxIter=200, stopCond=1e-7, initial=None, result="vecmat", retlogli=False):
    
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
