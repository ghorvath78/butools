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

def MMAPPH1NPPR(D, sigma, S, *argv):
    """
    Returns various performane measures of a continuous time 
    MMAP[K]/PH[K]/1 non-preemptive priority queue, see [1]_.
    
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
    
    # step 1. solution of the workload process of the joint queue
    # ===========================================================
    sM = np.sum(M)
    Qwmm = ml.matrix(D0)
    Qwpm = ml.zeros((N*sM, N))
    Qwmp = ml.zeros((N, N*sM))
    Qwpp = ml.zeros((N*sM, N*sM)) 
    kix = 0
    for i in range(K):
        Qwmp[:,kix:kix+N*M[i]] = np.kron(D[i+1], sigma[i])
        Qwpm[kix:kix+N*M[i],:] = np.kron(I,s[i])
        Qwpp[kix:kix+N*M[i],:][:,kix:kix+N*M[i]] = np.kron(I,S[i])
        kix += N*M[i]

    # calculate fundamental matrices
    Psiw, Kw, Uw = FluidFundamentalMatrices (Qwpp, Qwpm, Qwmp, Qwmm, 'PKU', precision)
    
    # calculate boundary vector
    Ua = ml.ones((N,1)) + 2*np.sum(Qwmp*(-Kw).I,1)
    pm = Linsolve (ml.hstack((Uw,Ua)).T, ml.hstack((ml.zeros((1,N)),ml.ones((1,1)))).T).T

    ro =  ((1.0-np.sum(pm))/2.0)/(np.sum(pm)+(1.0-np.sum(pm))/2.0) # calc idle time with weight=1, and the busy time with weight=1/2
    kappa = pm/np.sum(pm)
    
    pi = CTMCSolve (sD)
    lambd = []
    for i in range(K):
        lambd.append(np.sum(pi*D[i+1]))

    Psiw = []
    Qwmp = []
    Qwzp = []
    Qwpp = []
    Qwmz = []
    Qwpz = []
    Qwzz = []
    Qwmm = []
    Qwpm = []
    Qwzm = []
    for k in range(K):
        # step 2. construct a workload process for classes k...K
        # ======================================================
        Mlo = np.sum(M[:k])
        Mhi = np.sum(M[k:])

        Qkwpp = ml.zeros((N*Mlo*Mhi+N*Mhi, N*Mlo*Mhi+N*Mhi))
        Qkwpz = ml.zeros((N*Mlo*Mhi+N*Mhi, N*Mlo)) 
        Qkwpm = ml.zeros((N*Mlo*Mhi+N*Mhi, N))
        Qkwmz = ml.zeros((N, N*Mlo))
        Qkwmp = ml.zeros((N, N*Mlo*Mhi+N*Mhi))
        Dlo = ml.matrix(D0)
        for i in range(k):
            Dlo = Dlo + D[i+1]
        Qkwmm = Dlo
        Qkwzp = ml.zeros((N*Mlo, N*Mlo*Mhi+N*Mhi))
        Qkwzm = ml.zeros((N*Mlo, N))
        Qkwzz = ml.zeros((N*Mlo, N*Mlo))
        kix = 0
        for i in range(k,K):
            kix2 = 0
            for j in range(k):
                bs = N*M[j]*M[i]
                bs2 = N*M[j]
                Qkwpp[kix:kix+bs,kix:kix+bs] = np.kron(I,np.kron(ml.eye(M[j]),S[i]))
                Qkwpz[kix:kix+bs,kix2:kix2+bs2] = np.kron(I,np.kron(ml.eye(M[j]),s[i]))
                Qkwzp[kix2:kix2+bs2,kix:kix+bs] = np.kron(D[i+1],np.kron(ml.eye(M[j]), sigma[i]))
                kix += bs
                kix2 += bs2
        for i in range(k,K):
            bs = N*M[i]
            Qkwpp[kix:kix+bs,:][:,kix:kix+bs] = np.kron(I,S[i])
            Qkwpm[kix:kix+bs,:] = np.kron(I,s[i])
            Qkwmp[:,kix:kix+bs] = np.kron(D[i+1],sigma[i])
            kix += bs
        kix = 0
        for j in range(k):
            bs = N*M[j]
            Qkwzz[kix:kix+bs,kix:kix+bs] = np.kron(Dlo, ml.eye(M[j])) + np.kron(I, S[j])
            Qkwzm[kix:kix+bs,:] = np.kron(I, s[j])
            kix += bs

        if Qkwzz.shape[0]>0:
            Psikw = FluidFundamentalMatrices (Qkwpp+Qkwpz*(-Qkwzz).I*Qkwzp, Qkwpm+Qkwpz*(-Qkwzz).I*Qkwzm, Qkwmp, Qkwmm, 'P', precision)
        else:
            Psikw = FluidFundamentalMatrices (Qkwpp, Qkwpm, Qkwmp, Qkwmm, 'P', precision)
        Psiw.append(Psikw)
        
        Qwzp.append(Qkwzp)
        Qwmp.append(Qkwmp)
        Qwpp.append(Qkwpp)
        Qwmz.append(Qkwmz)
        Qwpz.append(Qkwpz)
        Qwzz.append(Qkwzz)
        Qwmm.append(Qkwmm)
        Qwpm.append(Qkwpm)
        Qwzm.append(Qkwzm)
    
    # step 3. calculate Phi vectors
    # =============================
    lambdaS = sum(lambd)
    phi = [(1-ro)*kappa*(-D0) / lambdaS]
    q0 = [[]]
    qL = [[]]
    for k in range(K-1):
        sDk = ml.matrix(D0)
        for j in range(k+1):
            sDk = sDk + D[j+1]
        # pk
        pk = sum(lambd[:k+1])/lambdaS - (1-ro)*kappa*np.sum(sDk,1)/lambdaS
        # A^(k,1)
        Qwzpk = Qwzp[k+1]
        vix = 0
        Ak = []
        for ii in range(k+1):
            bs = N*M[ii]
            V1 = Qwzpk[vix:vix+bs,:]
            Ak.append (np.kron(I,sigma[ii]) * (-np.kron(sDk,ml.eye(M[ii]))-np.kron(I,S[ii])).I * (np.kron(I,s[ii]) + V1*Psiw[k+1]))
            vix += bs
        # B^k
        Qwmpk = Qwmp[k+1]
        Bk = Qwmpk * Psiw[k+1]
        ztag = phi[0]*((-D0).I*D[k+1]*Ak[k] - Ak[0] + (-D0).I*Bk)
        for i in range(k):
            ztag += phi[i+1]*(Ak[i]-Ak[i+1]) + phi[0]*(-D0).I*D[i+1]*Ak[i]
        Mx = ml.eye(Ak[k].shape[0])-Ak[k]
        Mx[:,0] = ml.ones((N,1))
        phi.append(ml.hstack((pk, ztag[:,1:]))*Mx.I)  # phi(k) = Psi^(k)_k * p(k). Psi^(k)_i = phi(i) / p(k)

        q0.append(phi[0]*(-D0).I)
        qLii = []
        for ii in range(k+1):
            qLii.append((phi[ii+1] - phi[ii] + phi[0]*(-D0).I*D[ii+1]) * np.kron(I,sigma[ii]) * (-np.kron(sDk,ml.eye(M[ii]))-np.kron(I,S[ii])).I)
        qL.append(ml.hstack(qLii))
    
    
    # step 4. calculate performance measures
    # ======================================
    Ret = []
    for k in classes:

        sD0k = ml.matrix(D0)
        for i in range(k):
            sD0k +=  D[i+1]     
       
        if k<K-1:
            # step 4.1 calculate distribution of the workload process right 
            # before the arrivals of class k jobs
            # ============================================================
            if Qwzz[k].shape[0]>0:
                Kw = Qwpp[k]+Qwpz[k]*(-Qwzz[k]).I*Qwzp[k] + Psiw[k]*Qwmp[k]
            else:
                Kw = Qwpp[k] + Psiw[k]*Qwmp[k]
            BM = ml.zeros((0,0))
            CM = ml.zeros((0,N))
            DM = ml.zeros((0,0))
            for i in range(k):
                BM = la.block_diag(BM,np.kron(I,S[i]))
                CM = ml.vstack((CM, np.kron(I,s[i])))
                DM = la.block_diag(DM,np.kron(D[k+1],ml.eye(M[i])))
            if k>0:
                Kwu = ml.vstack((ml.hstack((Kw, (Qwpz[k]+Psiw[k]*Qwmz[k])*(-Qwzz[k]).I*DM)), ml.hstack((ml.zeros((BM.shape[0],Kw.shape[1])), BM))))
                Bwu = ml.vstack((Psiw[k]*D[k+1], CM))
                iniw = ml.hstack((q0[k]*Qwmp[k]+qL[k]*Qwzp[k], qL[k]*DM))
                pwu = q0[k]*D[k+1]
            else:
                Kwu = Kw
                Bwu = Psiw[k]*D[k+1]
                iniw = pm*Qwmp[k]
                pwu = pm*D[k+1]

            norm = np.sum(pwu) + np.sum(iniw*(-Kwu).I*Bwu)
            pwu = pwu / norm
            iniw = iniw / norm

            # step 4.2 create the fluid model whose first passage time equals the
            # WAITING time of the low prioroity customers
            # ==================================================================
            KN = Kwu.shape[0]
            Qspp = ml.zeros((KN+N*np.sum(M[k+1:]), KN+N*np.sum(M[k+1:])))
            Qspm = ml.zeros((KN+N*np.sum(M[k+1:]), N))
            Qsmp = ml.zeros((N, KN+N*np.sum(M[k+1:])))
            Qsmm = sD0k + D[k+1]
            kix = 0
            for i in range(k+1,K):
                bs = N*M[i]
                Qspp[KN+kix:KN+kix+bs,:][:,KN+kix:KN+kix+bs] = np.kron(I,S[i])
                Qspm[KN+kix:KN+kix+bs,:] = np.kron(I,s[i])
                Qsmp[:,KN+kix:KN+kix+bs] = np.kron(D[i+1],sigma[i])
                kix += bs

            Qspp[:KN,:][:,:KN] = Kwu
            Qspm[:KN,:] = Bwu
            inis = ml.hstack((iniw, ml.zeros((1,N*np.sum(M[k+1:])))))

            # calculate fundamental matrix
            Psis = FluidFundamentalMatrices (Qspp, Qspm, Qsmp, Qsmm, 'P', precision)

            # step 4.3. calculate the performance measures
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
                    # calculate waiting time moments
                    Pn = [Psis]
                    wtMoms = []
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
                        wtMoms.append(np.sum(inis*P*(-1)**n) / 2**n)
                    # calculate RESPONSE time moments
                    Pnr = [np.sum(inis*Pn[0])*sigma[k]]
                    rtMoms = []
                    for n in range(1,numOfSTMoms+1):
                        P =  n*Pnr[n-1]*(-S[k]).I + (-1)**n*np.sum(inis*Pn[n])*sigma[k] / 2**n
                        Pnr.append(P)
                        rtMoms.append(np.sum(P)+np.sum(pwu)*math.factorial(n)*np.sum(sigma[k]*(-S[k]).I**n))
                    Ret.append(rtMoms)
                    argIx += 1
                elif type(argv[argIx]) is str and argv[argIx]=="stDistr":
                    # DISTRIBUTION OF THE SOJOURN TIME
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    stCdfPoints = argv[argIx+1]
                    res = []
                    for t in stCdfPoints:
                        L = erlMaxOrder
                        lambdae = L/t/2
                        Psie = FluidFundamentalMatrices (Qspp-lambdae*ml.eye(Qspp.shape[0]), Qspm, Qsmp, Qsmm-lambdae*ml.eye(Qsmm.shape[0]), 'P', precision)
                        Pn = [Psie]
                        pr = (np.sum(pwu) + np.sum(inis*Psie)) * (1-np.sum(sigma[k]*(ml.eye(S[k].shape[0])-S[k]/2/lambdae).I**L))
                        for n in range(1,L):
                            A = Qspp + Psie*Qsmp - lambdae*ml.eye(Qspp.shape[0])
                            B = Qsmm + Qsmp*Psie - lambdae*ml.eye(Qsmm.shape[0])
                            C = 2*lambdae*Pn[n-1]
                            for i in range(1,n):
                                C += Pn[i]*Qsmp*Pn[n-i]
                            P = la.solve_sylvester(A, B, -C)
                            Pn.append(P)
                            pr += np.sum(inis*P) * (1-np.sum(sigma[k]*(np.eye(S[k].shape[0])-S[k]/2/lambdae).I**(L-n)))
                        res.append(pr)
                    Ret.append(np.array(res))
                    argIx += 1
                elif type(argv[argIx]) is str and (argv[argIx]=="ncMoms" or argv[argIx]=="ncDistr"):
                    W = (-np.kron(sD-D[k+1],ml.eye(M[k]))-np.kron(I,S[k])).I*np.kron(D[k+1],ml.eye(M[k]))
                    iW = (ml.eye(W.shape[0])-W).I
                    w = np.kron(ml.eye(N),sigma[k])
                    omega = (-np.kron(sD-D[k+1],ml.eye(M[k]))-np.kron(I,S[k])).I*np.kron(I,s[k])
                    if argv[argIx]=="ncMoms":
                        # MOMENTS OF THE NUMBER OF JOBS
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        numOfQLMoms = argv[argIx+1]
                        # first calculate it at departure instants
                        Psii = [Psis]
                        QLDPn = [inis*Psii[0]*w*iW]
                        for n in range(1,numOfQLMoms+1):
                            A = Qspp + Psis*Qsmp
                            B = Qsmm + Qsmp*Psis
                            C = n*Psii[n-1]*D[k+1]
                            bino = 1
                            for i in range(1,n):
                                bino = bino * (n-i+1) / i
                                C = C + bino * Psii[i]*Qsmp*Psii[n-i]
                            P = la.solve_sylvester(A, B, -C)
                            Psii.append(P)
                            QLDPn.append(n*QLDPn[n-1]*iW*W + inis*P*w*iW)
                        for n in range(numOfQLMoms+1):
                            QLDPn[n] = (QLDPn[n] + pwu*w*iW**(n+1)*W**n)*omega
                        # now calculate it at random time instance
                        QLPn = [pi]
                        qlMoms = []
                        iTerm = (ml.ones((N,1))*pi - sD).I
                        for n in range(1,numOfQLMoms+1):
                            sumP = np.sum(QLDPn[n]) + n*np.sum((QLDPn[n-1] - QLPn[n-1]*D[k+1]/lambd[k])*iTerm*D[k+1])
                            P = sumP*pi + n*(QLPn[n-1]*D[k+1] - QLDPn[n-1]*lambd[k])*iTerm
                            QLPn.append(P)
                            qlMoms.append(np.sum(P))
                        qlMoms = MomsFromFactorialMoms(qlMoms)
                        Ret.append(qlMoms)
                        argIx += 1
                    elif argv[argIx]=="ncDistr":
                        # DISTRIBUTION OF THE NUMBER OF JOBS
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        numOfQLProbs = argv[argIx+1]
                        Psid = FluidFundamentalMatrices (Qspp, Qspm, Qsmp, sD0k, 'P', precision)
                        Pn = [Psid]
                        XDn = inis*Psid*w
                        dqlProbs = (XDn+pwu*w)*omega
                        for n in range(1,numOfQLProbs):
                            A = Qspp + Psid*Qsmp
                            B = sD0k + Qsmp*Psid
                            C = Pn[n-1]*D[k+1]
                            for i in range(1,n):
                                C += Pn[i]*Qsmp*Pn[n-i]
                            P = la.solve_sylvester(A, B, -C)
                            Pn.append(P)
                            XDn = XDn*W + inis*P*w
                            dqlProbs = ml.vstack((dqlProbs, (XDn+pwu*w*W**n)*omega))
                        # now calculate it at random time instance
                        iTerm = -(sD-D[k+1]).I
                        qlProbs = lambd[k]*dqlProbs[0,:]*iTerm
                        for n in range(1,numOfQLProbs):
                            P = (qlProbs[n-1,:]*D[k+1]+lambd[k]*(dqlProbs[n,:]-dqlProbs[n-1,:]))*iTerm
                            qlProbs = ml.vstack((qlProbs, P))
                        qlProbs = np.sum(qlProbs,1).A.flatten()
                        Ret.append(qlProbs)
                        argIx += 1
                else:
                    raise Exception("MMAPPH1NPPR: Unknown parameter "+str(argv[argIx]))
                argIx += 1
        elif k==K-1:
            # step 3. calculate the performance measures
            # ==========================================   
            argIx = 0
            while argIx<len(argv):
                if argIx in eaten:
                    argIx += 1
                    continue
                elif type(argv[argIx]) is str and (argv[argIx]=="stMoms" or argv[argIx]=="stDistr"):
                    Kw = Qwpp[k]+Qwpz[k]*(-Qwzz[k]).I*Qwzp[k] + Psiw[k]*Qwmp[k]
                    AM = ml.zeros((0,0))
                    BM = ml.zeros((0,0))
                    CM = ml.zeros((0,1))
                    DM = ml.zeros((0,0))
                    for i in range(k):
                        AM = la.block_diag(AM,np.kron(ml.ones((N,1)),np.kron(ml.eye(M[i]),s[k])))
                        BM = la.block_diag(BM,S[i])
                        CM = ml.vstack((CM, s[i]))
                        DM = la.block_diag(DM,np.kron(D[k+1],ml.eye(M[i])))                        
                    Z = ml.vstack((ml.hstack((Kw, ml.vstack((AM,ml.zeros((N*M[k],AM.shape[1])))))), ml.hstack((ml.zeros((BM.shape[0],Kw.shape[1])), BM))))
                    z = ml.vstack((ml.zeros((AM.shape[0],1)), np.kron(ml.ones((N,1)),s[k]), CM))
                    iniw = ml.hstack((q0[k]*Qwmp[k]+qL[k]*Qwzp[k], ml.zeros((1,BM.shape[0]))))
                    zeta = iniw/np.sum(iniw*(-Z).I*z)
                    if argv[argIx]=="stMoms":
                        # MOMENTS OF THE SOJOURN TIME
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        numOfSTMoms = argv[argIx+1]
                        rtMoms = []
                        for i in range(1,numOfSTMoms+1):
                            rtMoms.append(np.sum(math.factorial(i)*zeta*(-Z).I**(i+1)*z))
                        Ret.append(rtMoms)
                        argIx += 1
                    if argv[argIx]=="stDistr":
                        # DISTRIBUTION OF THE SOJOURN TIME
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        stCdfPoints = argv[argIx+1]
                        rtDistr = []
                        for t in stCdfPoints:
                            rtDistr.append (np.sum(zeta*(-Z).I*(ml.eye(Z.shape[0])-la.expm(Z*t))*z))
                        Ret.append(np.array(rtDistr))
                        argIx += 1
                elif type(argv[argIx]) is str and (argv[argIx]=="ncMoms" or argv[argIx]=="ncDistr"):
                    L = ml.zeros((N*np.sum(M),N*np.sum(M)))
                    B = ml.zeros((N*np.sum(M),N*np.sum(M)))
                    F = ml.zeros((N*np.sum(M),N*np.sum(M)))
                    kix = 0
                    for i in range(K):
                        bs = N*M[i]
                        F[kix:kix+bs,:][:,kix:kix+bs] = np.kron(D[k+1],ml.eye(M[i]))
                        L[kix:kix+bs,:][:,kix:kix+bs] = np.kron(sD0k,ml.eye(M[i])) + np.kron(I,S[i])
                        if i<K-1:
                            L[kix:kix+bs,:][:,N*np.sum(M[:k]):] = np.kron(I,s[i]*sigma[k])
                        else:
                            B[kix:kix+bs,:][:,N*np.sum(M[:k]):] = np.kron(I,s[i]*sigma[k])
                        kix += bs
                    R = QBDFundamentalMatrices (B, L, F, 'R', precision)
                    p0 = ml.hstack((qL[k], q0[k]*np.kron(I,sigma[k])))
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
                    raise Exception("MMAPPH1NPPR: Unknown parameter "+str(argv[argIx]))
                argIx += 1

    if len(Ret)==1:
        return Ret[0]
    else:
        return Ret
