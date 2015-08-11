# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 08:24:03 2013

@author: Gabor Horvath
"""
import numpy as np
import numpy.matlib as ml
import math
from butools.mc import DTMCSolve

#def LikelihoodFromTrace (trace, A, B, prec=1e-14):
#
#    if A.shape[0]==1:
#        # uniformization based solution
#        tr = np.sort(trace)
#        trdiff = np.empty(len(tr))
#        trdiff[0] = tr[0]
#        trdiff[1:] = tr[1:]-tr[:-1]
#        lambd = np.max(np.abs(np.diag(B)))
#        loglambda = math.log(lambd)
#        P = B/lambd + ml.eye(B.shape[0])
#        a = np.sum(-B,1)
#        eps = max(prec, 10**(math.log(prec)/math.log(10.0) + math.log(lambd)/math.log(10.0)))
#        coeffv = ml.matrix(A)
#        logli = 0.0
#        for tk in trdiff:
#            if tk>0:
#                logtr = math.log(tk)
#                st = coeffv
#                lpoi = -lambd*tk
#                poi =  math.exp(lpoi)
#                spoi = poi
#                coeffv = st * poi
#                i = 1
#                while spoi<1.0-eps:
#                    lpoi += loglambda + logtr - math.log(i)
#                    poi =  math.exp(lpoi)
#                    spoi += poi
#                    st *= P
#                    coeffv += st*poi
#                    i += 1
#            logli += math.log(coeffv*a)
#        return logli/len(tr)
##        # naive solution -- extra slow!!!
##        l = 0
##        for t in trace:
##            l = l + math.log (np.sum(A*la.expm(B*t)*(-B)))
##        return l / len(trace)
#    else:
#        lambd = np.max(np.abs(np.diag(A)))
#        loglambda = math.log(lambd)       
#        P = A/lambd + ml.eye(A.shape[0])
#        eps = max(prec, 10**(math.log(prec)/math.log(10.0) + math.log(lambd)/math.log(10.0)))
#        coeffv = DTMCSolve(-A.I*B)
#        scale = 0
#        ix = 0
#        for tk in trace:
#            logtr = math.log(tk)
#            st = coeffv
#            lpoi = -lambd*tk
#            poi =  math.exp(lpoi)
#            spoi = poi
#            coeffv = st * poi
#            i = 1
#            while spoi<1.0-eps:
#                lpoi += loglambda + logtr - math.log(i)
#                poi =  math.exp(lpoi)
#                spoi += poi
#                st *= P
#                coeffv += st*poi
#                i += 1
#            coeffv *= B
#            csum = np.sum(coeffv)
#            while csum>1:
#                csum /= 2.0
#                coeffv /= 2.0
#                scale += 1
#            while csum<0.1:
#                csum *= 2.0
#                coeffv *= 2.0
#                scale -= 1
#            ix += 1
#        return (math.log(np.sum(coeffv)) + scale*math.log(2))/len(trace)

def LikelihoodFromTrace (trace, A, B, prec=1e-14):
    """
    Evaluates the log-likelihood of a trace with the given PH
    distribution or MAP. The result is divided by the length
    of the trace.
    
    If X is a row vector, than (X,Y) is interpreted as a PH
    distribution, otherwise (X,Y) is considered to be a MAP.
    
    Parameters
    ----------
    trace : column vector, length K
        The samples of the trace
    X : matrix, shape (1,M) or (M,M)
        If X is a row vector, it is the initial probability
        vector of the PH distribution. If X is a square
        matrix, it is interpreted as the D0 matrix of a MAP
    Y : matrix, (M,M)
        If X is a row vector, Y is the transient generator
        of the PH distribution. If X is a square matrix, Y
        is interpreted as the D1 matrix of a MAP
    prec : double, optional
        Numerical precision used by the randomization. The
        default value is 1e-14.
    
    Returns
    -------
    logli : double
        The log likelihood divided by the size of the trace
        
    Notes
    -----
    The procedure is much faster with PH distributions.
    """

    if A.shape[0]==1:
        # uniformization based solution
        tr = np.sort(trace)
        lambd = np.max(np.abs(np.diag(B)))
        loglambda = math.log(lambd)
        P = B/lambd + ml.eye(B.shape[0])
        a = np.sum(-B,1)
        eps = max(prec, 10**(math.log(prec)/math.log(10.0) + math.log(lambd)/math.log(10.0)))
        lpoi = -lambd*tr
        logtr = np.log(tr)
        poi =  np.exp(lpoi)
        spoi = np.array(poi)
        fx = poi*(A*a)[0,0]
        k = 1
        first = 0
        coeffv = ml.matrix(A)
        maxIter = 10000
        while k<maxIter:
            coeffv = coeffv * P
            lpoi[first:] += loglambda + logtr[first:] - math.log(k)
            poi[first:] = np.exp(lpoi[first:])
            spoi[first:] += poi[first:]
            fx[first:] += poi[first:] * (coeffv*a)[0,0]
            k += 1
            nfirst = (spoi[first:]<1-eps).nonzero()[0]
            if len(nfirst)==0:
                break
            first += nfirst[0]
        return np.sum(np.log(fx))/len(logtr)
    else:
        D0 = A
        D1 = B
        N = D0.shape[0]
        L = len(trace)
        
        # first we calculate matrix e^(D0*x(i))*D1 for each sample
        ix = np.argsort(trace)
        tr = trace[ix]
        lambd = np.max(np.abs(np.diag(D0)))
        loglambda = math.log(lambd)
        P = D0/lambd + ml.eye(N)
        eps = max(prec, 10**(math.log(prec)/math.log(10.0) + math.log(lambd)/math.log(10.0)))
        lpoi = -lambd*tr;
        logtr = np.log(tr)
        poi =  np.exp(lpoi)
        spoi = np.array(poi)
        coeffv = ml.matrix(D1)
        fx = np.kron(poi,coeffv)
        k = 1
        first = 0
        maxIter = 10000
        while k<maxIter:
            coeffv = P * coeffv
            lpoi[first:] += loglambda + logtr[first:] - math.log(k)
            poi[first:] = np.exp(lpoi[first:])
            spoi[first:] += poi[first:]           
            fx[:,first*N:] += np.kron(poi[first:],coeffv)
            k += 1 
            nfirst = (spoi[first:]<1-eps).nonzero()[0]
            if len(nfirst)==0:
                break
            first += nfirst[0]
        alpha = DTMCSolve ((-D0).I*D1)
        l = np.array(alpha)
        sc = 0
        ixrev = np.argsort(ix)
        for i in range(L):
            l = l.dot(fx[:,ixrev[i]*N:(ixrev[i]+1)*N])
            if i % 10 ==0:
                # sometimes we need to rescale the results to avoid "nan"s
                scale = math.ceil(math.log2(np.sum(l)))
                if scale>1:
                    l /= 2**scale
                    sc += scale
                if scale<-10:
                    scale += 10
                    l /= 2**scale
                    sc += scale
        return (math.log(np.sum(l))+sc*math.log(2)) / len(logtr)
