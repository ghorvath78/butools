# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 15:23:28 2013

@author: gabor
"""
import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
import numpy.matlib as ml
import math

def TransformToMonocyclic (A, maxSize=100, precision=1e-14):
    """
    Transforms an arbitrary matrix to a Markovian monocyclic 
    matrix (see [1]_).
    
    Parameters
    ----------
    A : matrix, shape (N,N)
        Matrix parameter of the initial representation
    maxSize : int, optional
        The maximal order of the resulting Markovian 
        representation. The default value is 100
    precision : double, optional
        Matrix entries smaller than the precision are 
        considered to be zeros. The default value is 1e-14
        
    Returns
    -------
    B : matrix, shape (M,M)
        Transient generator matrix of the Markovian monocyclic
        representation. Note that M>N if there are complex 
        eigenvalues.
    
    Notes
    -----    
    Raises an exception if no Markovian monocyclic generator 
    has been found up to the given size.
    
    References
    ----------
    .. [1]  Mocanu, S., Commault, C.: "Sparse representations 
            of phase-type distributions," Stoch. Models 15, 
            759-778 (1999)
    """
        
    def eigvalc(A):
        tol = math.sqrt(np.finfo(float).eps)
        lmb = np.sort(la.eigvals(A))
        lmb = np.around(lmb/tol) * tol
        evalues = np.unique(lmb)
        n = lmb.shape[0]
        d = evalues.shape[0]
        Ax = np.outer (np.ones((1,n)), evalues)
        Bx = np.outer (lmb, np.ones((1,d)))
        MATCH = np.abs(Ax-Bx) <= tol
        repeats = np.sum(MATCH,0)
        return (evalues, repeats)

    def generateFEBs (A, evalues, repeats, maxSize):

        # calculate the parameters of the febs
        febs = []
        i = 0
        size = 0
        while i<len(evalues):
            multip = repeats[i]
            feb = {}    
            evalimag = -abs(evalues[i].imag)
            if -evalimag < precision:
                feb = {"n": 1, "sigma": -evalues[i].real, "z": 0.0, "evals": np.array(evalues[i].real), "multip": multip}
                size += 1
            else:
                n = 3
                size += 3
                while evalimag / evalues[i].real >= 1.0 / math.tan(math.pi/n):
                    n += 1
                    size += 1
                    if size > maxSize:
                        raise Exception("The represetation is too large (>maxSize). No result returned.")
                sigma = -(2.0*evalues[i].real + evalimag*(1.0/math.tan(math.pi/n) - math.tan(math.pi/n))) / 2.0
                z = (-evalimag*(1.0/math.tan(math.pi/n) + math.tan(math.pi/n)) / (2.0*sigma))**n
                ev = []
                domi = evalues[i].real
                for k in range(1,n+1):
                    neval = complex(-(1.0-z**(1.0/n)*math.cos(2.0*(k-1)*math.pi/n))*sigma,  z**(1.0/n)*math.sin(2.0*(k-1)*math.pi/n)*sigma)
                    ev.append (neval)
                    if neval.real > domi:
                        domi = neval.real
                feb = {"n": n, "sigma": sigma, "z": z, "evals": np.array(ev), "multip": multip}
                
            febs.append(feb)
            if -evalimag < precision:
                i += 1
            else:
                i += 2
                
        # ordering according to the dominant eigenvalue (simple bubble sort)
        N = len(febs)
        for i in range(N):
            for j in range(N-i-1):
                if np.max(febs[j]["evals"]) < np.max(febs[j+1]["evals"]):
                    febs[j], febs[j+1] = febs[j+1], febs[j]
        return febs

    def febGenerator (lmb, z, n, multip):

        A = ml.zeros((multip*n,multip*n))
        for fi in range(multip):
            lv = np.ones(n)*lmb
            A[fi*n:(fi+1)*n, fi*n:(fi+1)*n] = -np.diag(lv) + np.diag(lv[0:n-1],1)
            A[(fi+1)*n-1,fi*n] += z * lmb
            if fi < multip-1:
                A[(fi+1)*n-1,(fi+1)*n] = (1-z)*lmb
        return A
        
    evalues, repeats = eigvalc (A)

    # assemble generator matrix of the hyper-feb 
    febs = generateFEBs (A, evalues, repeats, maxSize)
    fullN = 0
    for f in febs:
        fullN += f["n"]*f["multip"]
    
    G = ml.zeros((fullN,fullN))
    pos = 0
    for i in range(len(febs)):
        Ni = febs[i]["n"]*febs[i]["multip"]
        G[pos:pos+Ni,pos:pos+Ni] = febGenerator (febs[i]["sigma"], febs[i]["z"], febs[i]["n"], febs[i]["multip"])
        if i < len(febs)-1:
            G[pos+Ni-1, pos+Ni] = -np.sum(np.sum(G[pos+Ni-1,:],1),0)
        pos = pos + Ni

    return G

def TransformToAcyclic (A, maxSize=100, precision=1e-14):
    """
    Transforms an arbitrary matrix to a Markovian bi-diagonal 
    matrix.
    
    Parameters
    ----------
    A : matrix, shape (N,N)
        Matrix parameter of the initial representation
    maxSize : int, optional
        The maximal order of the resulting Markovian 
        representation. The default value is 100
    precision : double, optional
        Matrix entries smaller than the precision are 
        considered to be zeros. The default value is 1e-14
        
    Returns
    -------
    B : matrix, shape (N,N)
        Transient (bi-diagonal) generator matrix of the
        Markovian acyclic representation.
        
    Notes
    -----
    Calls the 'TransformToMonocyclic' procedure if all the 
    eigenvalues are real, otherwise it raises an error if no
    Markovian acyclic generator has been found up to the 
    given size.
    
    Raises an error if no Markovian acyclic generator 
    has been found up to the given size.
        References
    ----------
    .. [1]  Mocanu, S., Commault, C.: "Sparse representations 
            of phase-type distributions," Stoch. Models 15, 
            759-778 (1999)
    """

    eigs = la.eig (A)
    for ev in eigs:
        if np.any(np.abs(np.imag(ev))>=precision):
            raise Exception("TransformToAcyclic: Complex eigenvalue found, no acyclic representation exists.")
    
    return TransformToMonocyclic (A, maxSize, precision)

def ExtendToMarkovian (alpha, A, maxSize=100, precision=1e-14):
    """
    Assume we have an existing monocyclic (or acyclic) 
    representation of a matrix-exponential distribution 
    described by matrix A and vector alpha such that A is 
    Markovian but alpha is not. This procedure appends an 
    appropriate Erlang tail to the representation that makes 
    the result Markovian (both the generator matrix and the 
    initial vector parameter), while keeping the distribution 
    the same. In [1]_ it is proven that this is always 
    possible if the initial (alpha,A) defines a distribuion 
    (non-negative density).
    
    Parameters
    ----------
    alpha : vector, shape (1,M)
        The (non-Markovian) initial vector
    A : matrix, shape (M,M)            
        The (Markovian) transient generator.
    maxSize : int, optional
        The procedure stops if more than maxSize new phases 
        are required. The default value is 100
    precision : double, optional
        The initial vector is considered to be valid if the 
        smallest entry is greater than -precision. The
        default value is 1e-14
    
    Returns
    -------
    beta : vector, shape (1,N)
        The Markovian initial vector (N>=M)
    B : matrix, shape (N,N)
        The Markovian transient generator (N>=M).
    
    References
    ----------
    .. [1]  Mocanu, S., Commault, C.: "Sparse representations 
            of phase-type distributions," Stoch. Models 15, 
            759-778 (1999)
    """

    def addErlangTail (D, leng, mu):
        DN = D.shape[0]
        E = ml.zeros((DN+leng,DN+leng))
        E[0:DN,0:DN] = D
        E[DN-1,DN] = -np.sum(np.sum(D[DN-1,:],1),0)
        for ei in range(leng):
            E[DN+ei,DN+ei] = -mu
            if ei<leng-1:
                E[DN+ei,DN+ei+1] = mu
        return E

    def inivecWithTail (gamma, G, tailLength, mu):
        vlen = G.shape[0] + tailLength
        beta = ml.zeros((1, vlen))
        WG  = ml.eye(G.shape[0])+G/mu
        opv = np.copy(gamma)
        clv = -np.sum(G/mu,1)
        for k in range(vlen-1,G.shape[0]-1,-1):
            beta[0,k] = (opv * clv)[0,0]
            opv = opv * WG
        beta[0,0:G.shape[0]] = opv
        return beta

    # Step 1. find time shift t0 that makes the distribution positive
    # initial value of t0upper
    t0lower = 0.0
    t0upper = 1.0
    beta = alpha * sla.expm(A*t0upper)
    while np.min(beta) < -precision:
        t0upper *= 2.0
        beta = alpha * sla.expm(A*t0upper)
        
    # interval bisectioning
    while (t0upper - t0lower)/(t0upper + t0lower) > precision:
        t0 = (t0upper + t0lower) / 2.0
        beta = alpha * sla.expm(A*t0)
        if np.min(beta) < -precision:
            t0lower = t0
        else:
            t0upper = t0
    t0 = t0upper

    # find optimal length and rate parameters of the Erlang tail
    # points towards t0 sometimes give fewer states
    # thus we try to increase t0 gradually till the number of states
    # decreases
    
    increment = 1.1
    bestT0 = -1
    bestLupper = -1
    for i in range(100):
        # initial value of L0upper
        Llower = 1
        Lupper = 1
        beta = inivecWithTail (alpha, A, Lupper, Lupper/t0)
        while np.min(beta) < -precision and Lupper < maxSize-len(alpha):
            Lupper = Lupper * 2
            beta = inivecWithTail (alpha, A, Lupper, Lupper/t0)

        success =  np.min(beta)>=-precision
        if success:
            # interval bisectioning
            while Lupper - Llower > 1:
                L = int(round((Lupper + Llower) / 2))
                beta = inivecWithTail (alpha, A, L, L/t0)
                if np.min(beta) < -precision:
                    Llower = L
                else:
                    Lupper = L
        
        if success:
            if bestLupper>=0 and Lupper>bestLupper: # there was a successfull attempt before, and this one is worse
                break                               # stop and keep what we have
            else:                                   # otherwise go on and increment t0
                bestLupper = Lupper
                bestT0 = t0
                t0 *= increment
        else:
            if bestLupper>=0:
                break
            else:
                t0 *= increment

    t0 = bestT0
    Lupper = bestLupper
    
    if Lupper<0:
        raise Exception("ExtendToMarkovian: No positive representation found up to the given size!")        

    # final result
    beta = inivecWithTail (alpha, A, Lupper, Lupper/t0)
    B = addErlangTail (A, Lupper, Lupper/t0)
    return (beta, B)
    