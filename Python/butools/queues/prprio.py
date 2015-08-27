# -*- coding: utf-8 -*-
"""
Created on Tue May  5 11:31:18 2015

@author: gabor
"""
import numpy as np
import numpy.matlib as ml
import scipy.linalg as la
import butools
import math
from butools.moments import MomsFromFactorialMoms
from butools.map import CheckMMAPRepresentation
from butools.mc import CTMCSolve
from butools.ph import CheckPHRepresentation
from butools.mam import FluidFundamentalMatrices, QBDFundamentalMatrices
from butools.utils import Linsolve

def MMAPPH1PRPR(D, sigma, S, *argv):
    """
    Returns various performane measures of a MMAP[K]/PH[K]/1 
    preemptive resume priority queue, see [1]_.
    
    Parameters
    ----------
    D : list of matrices of shape (N,N), length (K+1)
        The D0...DK matrices of the arrival process.
        D1 corresponds to the lowest, DK to the highest priority.
    sigma : list of row vectors, length (K)
        The list containing the initial probability vectors of the service
        time distributions of the various customer types. The length of the
       vectors does not have to be the same.
    S : list of square matrices, length (K)
        The transient generators of the phase type distributions representing
        the service time of the jobs belonging to various types.
    further parameters : 
        The rest of the function parameters specify the options
        and the performance measures to be computed.
    
        The supported performance measures and options in this 
        function are:
    
        +----------------+--------------------+----------------------------------------+
        | Parameter name | Input parameters   | Output                                 |
        +================+====================+========================================+
        | "ncMoms"       | Number of moments  | The moments of the number of customers |
        +----------------+--------------------+----------------------------------------+
        | "ncDistr"      | Upper limit K      | The distribution of the number of      |
        |                |                    | customers from level 0 to level K-1    |
        +----------------+--------------------+----------------------------------------+
        | "stMoms"       | Number of moments  | The sojourn time moments               |
        +----------------+--------------------+----------------------------------------+
        | "stDistr"      | A vector of points | The sojourn time distribution at the   |
        |                |                    | requested points (cummulative, cdf)    |
        +----------------+--------------------+----------------------------------------+
        | "prec"         | The precision      | Numerical precision used as a stopping |
        |                |                    | condition when solving the Riccati and |
        |                |                    | the matrix-quadratic equations         |
        +----------------+--------------------+----------------------------------------+
        | "erlMaxOrder"  | Integer number     | The maximal Erlang order used in the   |
        |                |                    | erlangization procedure. The default   |
        |                |                    | value is 200.                          |
        +----------------+--------------------+----------------------------------------+
        | "classes"      | Vector of integers | Only the performance measures          |
        |                |                    | belonging to these classes are         |
        |                |                    | returned. If not given, all classes    |
        |                |                    | are analyzed.                          |
        +----------------+--------------------+----------------------------------------+
        
        (The quantities related to the number of customers in 
        the system include the customer in the server, and the 
        sojourn time related quantities include the service 
        times as well)
    
    Returns
    -------
    Ret : list of the performance measures
        Each entry of the list corresponds to a performance 
        measure requested. Each entry is a matrix, where the
        columns belong to the various job types.
        If there is just a single item, 
        then it is not put into a list.
    
    References
    ----------
    .. [1] G. Horvath, "Efficient analysis of the MMAP[K]/PH[K]/1
           priority queue", European Journal of Operational 
           Research, 246(1), 128-139, 2015.
    """
    
    K = len(D)-1

    # parse options
    eaten = []
    erlMaxOrder = 200;
    precision = 1e-14;
    classes = np.arange(0,K)
    for i in range(len(argv)):
        if argv[i]=="prec":
            precision = argv[i+1]
            eaten.append(i)
            eaten.append(i+1) 
        elif argv[i]=="erlMaxOrder":
            erlMaxOrder = argv[i+1]
            eaten.append(i)
            eaten.append(i+1) 
        elif argv[i]=="classes":
            classes = np.array(argv[i+1])-1
            eaten.append(i)
            eaten.append(i+1) 
    
    if butools.checkInput and not CheckMMAPRepresentation(D):
        raise Exception('MMAPPH1PRPR: The arrival process is not a valid MMAP representation!')
    
    if butools.checkInput:
        for k in range(K):
            if not CheckPHRepresentation(sigma[k],S[k]):
                raise Exception('MMAPPH1PRPR: the vector and matrix describing the service times is not a valid PH representation!')

    # some preparation
    D0 = D[0]
    N = D0.shape[0]
    I = ml.eye(N)
    sD = ml.zeros((N,N))
    for Di in D:
        sD += Di
    
    s = []
    M = np.empty(K)
    for i in range(K):
        s.append(np.sum(-S[i],1))
        M[i] = sigma[i].size
    
    Ret = []
    for k in classes:
       
        # step 1. solution of the workload process of the system
        # ======================================================
        sM = np.sum(M[k:K])
        Qwmm = ml.matrix(D0)
        for i in range(k):
            Qwmm += D[i+1]
            
        Qwpm = ml.zeros((N*sM, N))
        Qwmp = ml.zeros((N, N*sM))
        Qwpp = ml.zeros((N*sM, N*sM)) 
        kix = 0
        for i in range(k,K):
            Qwmp[:,kix:kix+N*M[i]] = np.kron(D[i+1], sigma[i])
            Qwpm[kix:kix+N*M[i],:] = np.kron(I,s[i])
            Qwpp[kix:kix+N*M[i],:][:,kix:kix+N*M[i]] = np.kron(I,S[i])
            kix += N*M[i]

        # calculate fundamental matrices
        Psiw, Kw, Uw = FluidFundamentalMatrices (Qwpp, Qwpm, Qwmp, Qwmm, 'PKU', precision)
        
        # calculate boundary vector
        Ua = ml.ones((N,1)) + 2*np.sum(Qwmp*(-Kw).I,1)
        pm = Linsolve (ml.hstack((Uw,Ua)).T, ml.hstack((ml.zeros((1,N)),ml.ones((1,1)))).T).T
      
        Bw = ml.zeros((N*sM, N))
        Bw[0:N*M[k],:] = np.kron(I,s[k])
        kappa = pm*Qwmp / np.sum(pm*Qwmp*(-Kw).I*Bw)
       
        if k<K-1:
            # step 2. construct fluid model for the remaining sojourn time process
            # ====================================================================
            # (for each class except the highest priority)
            Qsmm = ml.matrix(D0)
            for i in range(k+1):
                Qsmm += D[i+1]

            Np = Kw.shape[0]
            Qspm = ml.zeros((Np+N*np.sum(M[k+1:]), N))
            Qsmp = ml.zeros((N, Np+N*np.sum(M[k+1:])))
            Qspp = ml.zeros((Np+N*np.sum(M[k+1:]), Np+N*np.sum(M[k+1:])))
            Qspp[:Np,:Np] = Kw
            Qspm[:Np,:N] = Bw
            kix = Np
            for i in range(k+1,K):
                Qsmp[:,kix:kix+N*M[i]] = np.kron(D[i+1], sigma[i])
                Qspm[kix:kix+N*M[i],:] = np.kron(I,s[i])
                Qspp[kix:kix+N*M[i],kix:kix+N*M[i]] = np.kron(I,S[i])
                kix += N*M[i]

            inis = ml.hstack((kappa, ml.zeros((1,N*np.sum(M[k+1:])))))
            Psis = FluidFundamentalMatrices (Qspp, Qspm, Qsmp, Qsmm, 'P', precision)
            
            # step 3. calculate the performance measures
            # ==========================================   
            argIx = 0
            while argIx<len(argv):
                if argIx in eaten:
                    argIx += 1
                    continue
                elif type(argv[argIx]) is str and argv[argIx]=="stMoms":
                    # MOMENTS OF THE SOJOURN TIME
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    numOfSTMoms = argv[argIx+1]
                    Pn = [Psis]
                    rtMoms = []
                    for n in range(1,numOfSTMoms+1):
                        A = Qspp + Psis*Qsmp
                        B = Qsmm + Qsmp*Psis
                        C = -2*n*Pn[n-1]
                        bino = 1
                        for i in range(1,n):
                            bino = bino * (n-i+1) / i
                            C += bino * Pn[i]*Qsmp*Pn[n-i]
                        P = la.solve_sylvester(A, B, -C)
                        Pn.append(P)
                        rtMoms.append(np.sum(inis*P*(-1)**n) / 2**n)
                    Ret.append(rtMoms)
                    argIx += 1
                elif type(argv[argIx]) is str and argv[argIx]=="stDistr":
                    # DISTRIBUTION OF THE SOJOURN TIME
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    stCdfPoints = argv[argIx+1]
                    res = []
                    for t in stCdfPoints:
                        L = erlMaxOrder
                        lambd = L/t/2
                        Psie = FluidFundamentalMatrices (Qspp-lambd*ml.eye(Qspp.shape[0]), Qspm, Qsmp, Qsmm-lambd*ml.eye(Qsmm.shape[0]), 'P', precision)
                        Pn = [Psie]
                        pr = np.sum(inis*Psie)
                        for n in range(1,L):
                            A = Qspp + Psie*Qsmp - lambd*ml.eye(Qspp.shape[0])
                            B = Qsmm + Qsmp*Psie - lambd*ml.eye(Qsmm.shape[0])
                            C = 2*lambd*Pn[n-1]
                            for i in range(1,n):
                                C += Pn[i]*Qsmp*Pn[n-i]
                            P = la.solve_sylvester(A, B, -C)
                            Pn.append(P)
                            pr += np.sum(inis*P)
                        res.append(pr)
                    Ret.append(np.array(res))
                    argIx += 1
                elif type(argv[argIx]) is str and argv[argIx]=="ncMoms":
                    # MOMENTS OF THE NUMBER OF JOBS
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    numOfQLMoms = argv[argIx+1]
                    # first calculate it at departure instants
                    QLDPn = [Psis]
                    dqlMoms = []
                    for n in range(1,numOfQLMoms+1):
                        A = Qspp + Psis*Qsmp
                        B = Qsmm + Qsmp*Psis
                        C = n*QLDPn[n-1]*D[k+1]
                        bino = 1
                        for i in range(1,n):
                            bino = bino * (n-i+1) / i
                            C = C + bino * QLDPn[i]*Qsmp*QLDPn[n-i]
                        P = la.solve_sylvester(A, B, -C)
                        QLDPn.append(P)
                        dqlMoms.append(np.sum(inis*P))
                    dqlMoms = MomsFromFactorialMoms(dqlMoms)
                    # now calculate it at random time instance
                    pi = CTMCSolve (sD)
                    lambdak = np.sum(pi*D[k+1])
                    QLPn = [pi]
                    qlMoms = []
                    iTerm = (ml.ones((N,1))*pi - sD).I
                    for n in range(1,numOfQLMoms+1):
                        sumP = np.sum(inis*QLDPn[n]) + n*(inis*QLDPn[n-1] - QLPn[n-1]*D[k+1]/lambdak)*iTerm*np.sum(D[k+1],1)
                        P = sumP*pi + n*(QLPn[n-1]*D[k+1] - inis*QLDPn[n-1]*lambdak)*iTerm
                        QLPn.append(P)
                        qlMoms.append(np.sum(P))
                    qlMoms = MomsFromFactorialMoms(qlMoms)
                    Ret.append(qlMoms)
                    argIx += 1
                elif type(argv[argIx]) is str and argv[argIx]=="ncDistr":
                    # DISTRIBUTION OF THE NUMBER OF JOBS
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    numOfQLProbs = argv[argIx+1]
                    sDk = ml.matrix(D0)
                    for i in range(k):
                        sDk +=  D[i+1]
                    # first calculate it at departure instants
                    Psid = FluidFundamentalMatrices (Qspp, Qspm, Qsmp, sDk, 'P', precision)
                    Pn = [Psid]
                    dqlProbs = inis*Psid
                    for n in range(1,numOfQLProbs):
                        A = Qspp + Psid*Qsmp
                        B = sDk + Qsmp*Psid
                        C = Pn[n-1]*D[k+1]
                        for i in range(1,n):
                            C += Pn[i]*Qsmp*Pn[n-i]
                        P = la.solve_sylvester(A, B, -C)
                        Pn.append(P)
                        dqlProbs = ml.vstack((dqlProbs, inis*P))
                    # now calculate it at random time instance
                    pi = CTMCSolve (sD)
                    lambdak = np.sum(pi*D[k+1])
                    iTerm = -(sD-D[k+1]).I
                    qlProbs = lambdak*dqlProbs[0,:]*iTerm
                    for n in range(1,numOfQLProbs):
                        P = (qlProbs[n-1,:]*D[k+1]+lambdak*(dqlProbs[n,:]-dqlProbs[n-1,:]))*iTerm
                        qlProbs = ml.vstack((qlProbs, P))
                    qlProbs = np.sum(qlProbs,1).A.flatten()
                    Ret.append(qlProbs)
                    argIx += 1
                else:
                    raise Exception("MMAPPH1PRPR: Unknown parameter "+str(argv[argIx]))
                argIx += 1
        elif k==K-1:
            # step 3. calculate the performance measures
            # ==========================================   
            argIx = 0
            while argIx<len(argv):
                if argIx in eaten:
                    argIx += 1
                    continue
                elif type(argv[argIx]) is str and argv[argIx]=="stMoms":
                    # MOMENTS OF THE SOJOURN TIME
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    numOfSTMoms = argv[argIx+1]
                    rtMoms = []
                    for i in range(1,numOfSTMoms+1):
                        rtMoms.append(np.sum(math.factorial(i)*kappa*(-Kw).I**(i+1)*Bw))
                    Ret.append(rtMoms)
                    argIx += 1
                elif type(argv[argIx]) is str and argv[argIx]=="stDistr":
                    # DISTRIBUTION OF THE SOJOURN TIME
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    stCdfPoints = argv[argIx+1]
                    rtDistr = []
                    for t in stCdfPoints:
                        rtDistr.append (np.sum(kappa*(-Kw).I*(ml.eye(Kw.shape[0])-la.expm(Kw*t))*Bw))
                    Ret.append(np.array(rtDistr))
                    argIx += 1
                elif type(argv[argIx]) is str and (argv[argIx]=="ncMoms" or argv[argIx]=="ncDistr"):
                    L = np.kron(sD-D[k+1],ml.eye(M[k]))+np.kron(ml.eye(N),S[k])
                    B = np.kron(ml.eye(N),s[k]*sigma[k])
                    F = np.kron(D[k+1],ml.eye(M[k]))
                    L0 = np.kron(sD-D[k+1],ml.eye(M[k]))
                    R = QBDFundamentalMatrices (B, L, F, 'R', precision)
                    p0 = CTMCSolve(L0+R*B)
                    p0 = p0/np.sum(p0*(ml.eye(R.shape[0])-R).I)
                    if argv[argIx]=="ncMoms":
                        # MOMENTS OF THE NUMBER OF JOBS
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        numOfQLMoms = argv[argIx+1]
                        qlMoms = []
                        for i in range(1,numOfQLMoms+1):
                            qlMoms.append(np.sum(math.factorial(i)*p0*R**i*(ml.eye(R.shape[0])-R).I**(i+1)))
                        Ret.append(MomsFromFactorialMoms(qlMoms))
                    elif argv[argIx]=="ncDistr":
                        # DISTRIBUTION OF THE NUMBER OF JOBS
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        numOfQLProbs = argv[argIx+1]
                        qlProbs = [np.sum(p0)]
                        for i in range(1,numOfQLProbs):
                            qlProbs.append(np.sum(p0*R**i))
                        Ret.append(np.array(qlProbs))
                    argIx += 1
                else:
                    raise Exception("MMAPPH1PRPR: Unknown parameter "+str(argv[argIx]))
                argIx += 1

    if len(Ret)==1:
        return Ret[0]
    else:
        return Ret

