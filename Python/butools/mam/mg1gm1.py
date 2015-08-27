# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 08:33:12 2013

@author: gabor
"""

import numpy as np
import numpy.matlib as ml
import numpy.linalg as la
from scipy.fftpack import fft, ifft
import butools
from butools.mc import DTMCSolve, DRPSolve
from butools.utils import Diag
import math
import cmath

def GM1TypeCaudal(A):
    """
    Computes the Spectral Radius of R.

    [eta,v]=GIM1_CAUDAL(A) computes the dominant eigenvalue of the
    matrix R, the smallest nonnegative solution to 
    R= A0 + R A1 + R^2 A2 + ... + R^max Amax, if the GI/M/1
    type Markov chain characterized by A is recurrent.
    eta is the unique solution of PF(A0 + A1 z + A2 z^2 + ... 
    + Amax z^max) = z on (0,1), where PF denotes the Peron-Frobenius 
    eigenvalue. The right eigenvector v corresponding to the
    Peron-Frobenius eigenvalue of A(eta) is also returned.  

    Notes
    -----
    This procedure is the direct port of "SMCSolver: A MATLAB Toolbox for solving 
    Quasi-Birth-and-Death (QBD) type Markov chains". When using this procedure, 
    please refer to the paper Structured Markov chains solver: software tools by 
    Bini, Meini, Steffe and Van Houdt (SMCtools workshop).   
    """
    m = A.shape[0]
    dega = A.shape[1]//m-1
    eta_min = 0
    eta_max = 1
    eta = 0.5
    while eta_max - eta_min > 1e-14:
        temp = A[:,dega*m:]
        for i in range(dega-1,-1,-1):
            temp = temp * eta + A[:,i*m:(i+1)*m]
        new_eta = np.max(la.eigvals(temp))
        if new_eta > eta:
            eta_min = eta
        else:
            eta_max = eta
        eta = (eta_min+eta_max) / 2.0

    D, V = la.eig(temp)
    return (eta.real, np.real(V[:,np.argmax(D)]))

def MG1TypeDecay(A):
    """
    Computes the Decay Rate of the MG1 type MC [Falkenberg]

    [eta,uT]=MG1_Decay(A) computes the decay rate of a recurrent M/G/1 
    type Markov chain and the left eigenvector corresponding to the
    Peron-Frobenius eigenvalue of A(eta).
    Eta is the unique solution of PF(A0 + A1 z + A2 z^2 + ... + Amax z^max) = z on (1,RA), 
    where PF denotes the Peron-Frobenius eigenvalue. 

    Notes
    -----
    This procedure is the direct port of "SMCSolver: A MATLAB Toolbox for solving 
    Quasi-Birth-and-Death (QBD) type Markov chains". When using this procedure, 
    please refer to the paper Structured Markov chains solver: software tools by 
    Bini, Meini, Steffe and Van Houdt (SMCtools workshop).
    """
    m = A.shape[0]
    dega = A.shape[1]/m-1;
    eta = 1
    new_eta = 0
    while new_eta - eta < 0:
        eta += 1
        temp = A[:,dega*m:]
        for i in range(dega-1,-1,-1):
            temp = temp * eta + A[:,i*m:(i+1)*m]
        new_eta = np.max(la.eigvals(temp))

    eta_min = eta-1
    eta_max = eta
    eta = eta_min+0.5
    while eta_max - eta_min > 1e-14:
        temp = A[:,dega*m:]
        for i in range(dega-1,-1,-1):
            temp = temp * eta + A[:,i*m:(i+1)*m]
        new_eta = np.max(la.eigvals(temp))
        if new_eta < eta:
            eta_min = eta;
        else:
            eta_max = eta;
        eta = (eta_min+eta_max) / 2.0;

    D, V = la.eig(temp.T)
    return (eta.real, np.real(V[:,np.argmax(D)]).T)

def MG1TypeShifts (A, shiftType):
    """
    This function performs the shift operation for the MG1 case
    input: A = [A0 A1 A2 ... A_maxd]
    includes: one, tau and dbl

    Notes
    -----
    This procedure is the direct port of "SMCSolver: A MATLAB Toolbox for solving 
    Quasi-Birth-and-Death (QBD) type Markov chains". When using this procedure, 
    please refer to the paper Structured Markov chains solver: software tools by 
    Bini, Meini, Steffe and Van Houdt (SMCtools workshop).
    """
    m = A.shape[0]
    I = ml.eye(m)
    v = ml.zeros((m,1)) # default value if redundant
    tau = 1 # default value if redundant
    maxd = A.shape[1]//m-1
    sumA = ml.matrix(A[:,maxd*m:])
    beta = np.sum(sumA,1)
    # beta = (A_maxd)e + (A_maxd + A_maxd-1)e + ... + (Amaxd+...+A1)e
    for i in range(maxd-1,0,-1):
        sumA += A[:,i*m:(i+1)*m]
        beta += np.sum(sumA,1)
    sumA += A[:,:m]
    theta = DTMCSolve(sumA)
    drift = (theta * beta)[0,0]
    
    # if drift < 1 : positive recurrent
    hatA = ml.zeros((m,m*(maxd+1)))
    if drift < 1:
        if shiftType=="tau" or shiftType=="dbl":  # shift tau to infinity
            tau, uT = MG1TypeDecay(A)
            A[:,m:2*m] = A[:,m:2*m]-I
            uT = uT / np.sum(uT)
            rowhatA = ml.zeros((1,m*(maxd+1)))
            rowhatA[0,maxd*m:] = uT*A[:,maxd*m:]
            for i in range(maxd-1,-1,-1):
                rowhatA[0,i*m:(i+1)*m] = tau*rowhatA[0,(i+1)*m:(i+2)*m] + uT*A[:,i*m:(i+1)*m]
            hatA = A - ml.ones((m,1)) * rowhatA
        if shiftType=="dbl": # shift one to zero
            A = hatA
            # e is also the right eigenvector of hatA(1)
            # as the shift-tau does not influence G and Ge=e,
            # implying that hatA(1)e = e as G=hatA(G)
        if shiftType=="one": # shift one to zero
            A[:,m:2*m] = A[:,m:2*m] - I
        if shiftType=="one" or shiftType=="dbl": # shift one ot zero
            colhatA = ml.zeros((m,maxd+1))
            colhatA[:,0] = np.sum(A[:,0:m],1) # colhatA(:,1) = (A0)e
            for i in range(1,maxd+1):
                colhatA[:,i] = colhatA[:,i-1] + np.sum(A[:,i*m:(i+1)*m],1)
                # colhatA(:,i+1) = (A0+A1+...+Ai)e
            hatA = A - np.kron(colhatA, ml.ones((1,m))/m) # hatAi = Ai - (A0+A1+...+Ai)e*uT
    else:
        if shiftType=="one" or shiftType=="dbl": # shift one to infinity
            A[:,m:2*m] = A[:,m:2*m] - I
            rowhatA = ml.zeros((1,m*(maxd+1)))
            rowhatA[0,maxd*m:] = theta*A[:,maxd*m:]
            for i in range(maxd-1,-1,-1):
                rowhatA[0,i*m:(i+1)*m] = rowhatA[0,(i+1)*m:(i+2)*m] + theta*A[:,i*m:(i+1)*m] # rowhatAi = theta(Amaxd+...+Ai)
            hatA = A - ml.ones((m,1)) * rowhatA
        if shiftType=="dbl": # shift one to infinity
            A = hatA
            A[:,m:2*m] = A[:,m:2*m] + I
            # v is also the right eigenvector of hatA(tau)
            # as the shift-one does not influence G and Gv=tau*v,
            # implying that hatA(tau)v = tau*v as G=hatA(G)
        if shiftType=="tau" or shiftType=="dbl":  # shift tau to zero
            tau,v = GM1TypeCaudal(A)
            A[:,m:2*m] = A[:,m:2*m] - I
            v = v/np.sum(v)
            colhatA = ml.zeros((m,maxd+1))
            colhatA[:,0] = A[:,0:m] * v # colhatA(:,1) = (A0)v
            for i in range(1,maxd+1):
                colhatA[:,i] = colhatA[:,i-1]*tau**(-1) + A[:,i*m:(i+1)*m]*v
            hatA = A - np.kron(colhatA, ml.ones((1,m))) 
    hatA[:,m:2*m] = hatA[:,m:2*m] + I
    return (hatA,drift,tau,v)

def MG1FundamentalMatrix (A, precision=1e-14, maxNumIt=50, method="ShiftPWCR", maxNumRoot=2048, shiftType="one"):
    """
    Returns matrix G corresponding to the M/G/1 type Markov
    chain defined by matrices A.
    
    Matrix G is the minimal non-negative solution of the 
    following matrix equation:
    
    .. math::
        G = A_0 + A_1 G + A_2 G^2 + A_3 G^3 + \dots.
    
    The implementation is based on [1]_, please cite it if
    you use this method.
    
    Parameters
    ----------
    A : length(M) list of matrices of shape (N,N)
        Matrix blocks of the M/G/1 type generator from
        0 to M-1.
    precision : double, optional
        Matrix G is computed iteratively up to this
        precision. The default value is 1e-14
    maxNumIt : int, optional
        The maximal number of iterations. The default value
        is 50.
    method : {"CR", "RR", "NI", "FI", "IS"}, optional
        The method used to solve the matrix-quadratic
        equation (CR: cyclic reduction, RR: Ramaswami
        reduction, NI: Newton iteration, FI: functional
        iteration, IS: invariant subspace method). The 
        default is "CR".
    
    Returns
    -------
    G : matrix, shape (N,N)
        The G matrix of the M/G/1 type Markov chain.
        (G is stochastic.)
    
    References
    ----------
    .. [1] Bini, D. A., Meini, B., Steffé, S., Van Houdt,
           B. (2006, October). Structured Markov chains 
           solver: software tools. In Proceeding from the
           2006 workshop on Tools for solving structured 
           Markov chains (p. 14). ACM.
    """

    if not isinstance(A,np.ndarray):
        D = np.hstack(A)
    else:
        D = ml.matrix(A)
   
    if method=="ShiftPWCR":
        if butools.verbose:
            Dold = ml.matrix(D)
        D, drift, tau, v = MG1TypeShifts (D, shiftType)

    m = D.shape[0]
    M = D.shape[1]
    I = ml.eye(m)
    D = D.T
    D = np.vstack((D, ml.zeros(((2**(1+math.floor(math.log(M/m-1,2)))+1)*m-M,m))))

    # Step 0
    G = ml.zeros((m,m))
    Aeven = D[np.remainder(np.kron(np.arange(D.shape[0]/m),np.ones(m)),2)==0,:]
    Aodd  = D[np.remainder(np.kron(np.arange(D.shape[0]/m),np.ones(m)),2)==1,:]
    
    Ahatodd = np.vstack((Aeven[m:,:], D[-m:,:]))
    Ahateven = Aodd

    Rj = D[m:2*m,:]
    for i in range(3,M//m+1):
        Rj = Rj + D[(i-1)*m:i*m,:]
    Rj = la.inv(I-Rj)
    Rj = D[:m,:] * Rj
    
    numit = 0

    while numit < maxNumIt:
        numit+=1
        nj=Aodd.shape[0]//m
        if nj > 0:
            # Evaluate the 4 functions in the nj+1 roots using FFT
            # prepare for FFTs (such that they can be performed in 4 calls)
            
            temp1=np.reshape(Aodd[:nj*m,:].T,(m*m,nj),order='F').T
            temp2=np.reshape(Aeven[:nj*m,:].T,(m*m,nj),order='F').T
            temp3=np.reshape(Ahatodd[:nj*m,:].T,(m*m,nj),order='F').T
            temp4=np.reshape(Ahateven[:nj*m,:].T,(m*m,nj),order='F').T
            # FFTs           
            temp1=fft(temp1,nj,0)
            temp2=fft(temp2,nj,0)
            temp3=fft(temp3,nj,0)
            temp4=fft(temp4,nj,0)            
            # reform the 4*nj matrices
            temp1=np.reshape(temp1.T,(m,m*nj),order='F').T
            temp2=np.reshape(temp2.T,(m,m*nj),order='F').T
            temp3=np.reshape(temp3.T,(m,m*nj),order='F').T
            temp4=np.reshape(temp4.T,(m,m*nj),order='F').T
            
            # Next, we perform a point-wise evaluation of (6.20) - Thesis Meini
            Ahatnew = ml.empty((nj*m,m), dtype=complex)
            Anew = ml.empty((nj*m,m), dtype=complex)
            for cnt in range(1,nj+1):
                Ahatnew[(cnt-1)*m:cnt*m,:] = temp4[(cnt-1)*m:cnt*m,:] + temp2[(cnt-1)*m:cnt*m,:] * la.inv(I-temp1[(cnt-1)*m:cnt*m,:]) * temp3[(cnt-1)*m:cnt*m,:]
                Anew[(cnt-1)*m:cnt*m,:] = cmath.exp(-(cnt-1)*2.0j*math.pi/nj) * temp1[(cnt-1)*m:cnt*m,:] + temp2[(cnt-1)*m:cnt*m,:] * la.inv(I-temp1[(cnt-1)*m:cnt*m,:]) * temp2[(cnt-1)*m:cnt*m,:]
    
            # We now invert the FFTs to get Pz and Phatz   
            # prepare for IFFTs (in 2 calls)
            Ahatnew = np.reshape(Ahatnew[:nj*m,:].T,(m*m,nj),order='F').T
            Anew = np.reshape(Anew[:nj*m,:].T,(m*m,nj),order='F').T    
            # IFFTs
            Ahatnew = np.real(ifft(Ahatnew,nj,0))
            Anew = np.real(ifft(Anew,nj,0))
            # reform matrices Pi and Phati
            Ahatnew = np.reshape(Ahatnew.T,(m,m*nj),order='F').T
            Anew = np.reshape(Anew.T,(m,m*nj),order='F').T            
        else: # series Aeven, Aodd, Ahateven and Ahatodd are constant
            temp = Aeven * la.inv(I-Aodd)
            Ahatnew = Ahateven + temp*Ahatodd
            Anew = np.vstack((temp*Aeven, Aodd))
    
        nAnew = 0
        deg = Anew.shape[0]//m
        for i in range(deg//2,deg):
            nAnew = max(nAnew, la.norm(Anew[i*m:(i+1)*m,:],np.inf))

        nAhatnew = 0
        deghat = Ahatnew.shape[0]//m
        for i in range(deghat//2,deghat):
            nAhatnew = max(nAhatnew, la.norm(Ahatnew[i*m:(i+1)*m,:],np.inf))
        
        # c) the test
        while (nAnew > nj*precision or nAhatnew > nj*precision) and nj < maxNumRoot:
    
            nj *= 2
            stopv = min(nj, Aodd.shape[0]/m)
    
            # prepare for FFTs
            temp1=np.reshape(Aodd[:stopv*m,:].T,(m*m,stopv),order='F').T
            temp2=np.reshape(Aeven[:stopv*m,:].T,(m*m,stopv),order='F').T
            temp3=np.reshape(Ahatodd[:stopv*m,:].T,(m*m,stopv),order='F').T
            temp4=np.reshape(Ahateven[:stopv*m,:].T,(m*m,stopv),order='F').T
            # FFTs
            temp1=fft(temp1,nj,0)
            temp2=fft(temp2,nj,0)
            temp3=fft(temp3,nj,0)
            temp4=fft(temp4,nj,0)
            # reform the 4*(nj+1) matrices
            temp1=np.reshape(temp1.T,(m,m*nj),order='F').T
            temp2=np.reshape(temp2.T,(m,m*nj),order='F').T
            temp3=np.reshape(temp3.T,(m,m*nj),order='F').T
            temp4=np.reshape(temp4.T,(m,m*nj),order='F').T
    
            # Next, we perform a point-wise evaluation of (6.20) - Thesis Meini
            Ahatnew = ml.empty((nj*m,m), dtype=complex)
            Anew = ml.empty((nj*m,m), dtype=complex)
            for cnt in range(1,nj+1):
                Ahatnew[(cnt-1)*m:cnt*m,:] = temp4[(cnt-1)*m:cnt*m,:] + temp2[(cnt-1)*m:cnt*m,:] * la.inv(I-temp1[(cnt-1)*m:cnt*m,:]) * temp3[(cnt-1)*m:cnt*m,:]
                Anew[(cnt-1)*m:cnt*m,:] = cmath.exp(-(cnt-1)*2j*math.pi/nj) * temp1[(cnt-1)*m:cnt*m,:] + temp2[(cnt-1)*m:cnt*m,:] * la.inv(I-temp1[(cnt-1)*m:cnt*m,:]) * temp2[(cnt-1)*m:cnt*m,:]

            # We now invert the FFTs to get Pz and Phatz
            # prepare for IFFTs
            Ahatnew = np.reshape(Ahatnew[:nj*m,:].T,(m*m,nj),order='F').T
            Anew = np.reshape(Anew[:nj*m,:].T,(m*m,nj),order='F').T   
            # IFFTs
            Ahatnew = ml.matrix(np.real(ifft(Ahatnew,nj,0)))
            Anew = ml.matrix(np.real(ifft(Anew,nj,0)))
            # reform matrices Pi and Phati
            Ahatnew = np.reshape(Ahatnew.T,(m,m*nj),order='F').T
            Anew = np.reshape(Anew.T,(m,m*nj),order='F').T
    
            vec1 = ml.zeros((1,m))
            vec2 = ml.zeros((1,m))
            for i in range(1,Anew.shape[0]//m):
                vec1 += i*np.sum(Anew[i*m:(i+1)*m,:],0)
                vec2 += i*np.sum(Ahatnew[i*m:(i+1)*m,:],0)

            nAnew = 0
            deg = Anew.shape[0]//m
            for i in range(deg//2,deg):
                nAnew = max(nAnew, la.norm(Anew[i*m:(i+1)*m,:],np.inf))

            nAhatnew = 0
            deghat = Ahatnew.shape[0]//m
            for i in range(deghat//2,deghat):
                nAhatnew = max(nAhatnew, la.norm(Ahatnew[i*m:(i+1)*m,:],np.inf))
        if (nAnew > nj*precision or nAhatnew > nj*precision) and nj >= maxNumRoot:
            print("MaxNumRoot reached, accuracy might be affected!")
    
        if nj > 2:
            Anew = Anew[:m*nj/2,:]
            Ahatnew = Ahatnew[:m*nj/2,:]
        
        # compute Aodd, Aeven, ...
        Aeven = Anew[np.remainder(np.kron(np.arange(Anew.shape[0]/m),np.ones(m)),2)==0,:]
        Aodd = Anew[np.remainder(np.kron(np.arange(Anew.shape[0]/m),np.ones(m)),2)==1,:]
        
        Ahateven = Ahatnew[np.remainder(np.kron(np.arange(Ahatnew.shape[0]/m),np.ones(m)),2)==0,:]
        Ahatodd = Ahatnew[np.remainder(np.kron(np.arange(Ahatnew.shape[0]/m),np.ones(m)),2)==1,:]
        
        if butools.verbose==True:
            if method == "PWCR":
                print("The Point-wise evaluation of Iteration ", numit, " required ", nj, " roots")
            else:
                print("The Shifted PWCR evaluation of Iteration ", numit, " required ", nj, " roots")
        
        # test stopcriteria
        if method == "PWCR":
            Rnewj = Anew[m:2*m,:]
            for i in range (3, Anew.shape[0]/m+1):
                Rnewj = Rnewj + Anew[(i-1)*m:i*m,:]
            Rnewj = la.inv(I-Rnewj)
            Rnewj = Anew[:m,:]*Rnewj

            if np.max(np.abs(Rj-Rnewj)) < precision or np.max(np.sum(I-Anew[:m,:]*la.inv(I-Anew[m:2*m,:]),0)) < precision:
                G = Ahatnew[:m,:]
                for i in range (2,Ahatnew.shape[0]/m+1):
                    G = G + Rnewj * Ahatnew[(i-1)*m:i*m,:]
                G = D[:m,:]*la.inv(I-G)
                break
            
            Rj = Rnewj
            # second condition tests whether Ahatnew is degree 0 (numerically)
            if la.norm(Anew[:m,:m]) < precision or np.sum(Ahatnew[m:,:]) < precision or np.max(np.sum(I-D[:m,:]*la.inv(I-Ahatnew[:m,:]),0)) < precision:
                G = D[0:m,:] * la.inv(I-Ahatnew[:m,:])
                break
        else:
            Gold = G
            G = D[:m,:]*la.inv(I-Ahatnew[:m,:])
            if la.norm(G-Gold,np.inf) < precision or la.norm(Ahatnew[m:,:],np.inf) < precision:
                break
    
    if numit == maxNumIt and G==ml.zeros((m,m)):
        print("Maximum Number of Iterations reached!")
        G = D[:m,:] * la.inv(I-Ahatnew[:m,:])
    
    G=G.T
    
    if method=="ShiftPWCR":
        if shiftType=="one":
            G = G + (drift<1)*ml.ones((m,m))/m
        elif shiftType=="tau":
            G = G + (drift>1)*tau*v*ml.ones((1,m))
        elif shiftType=="dbl":
            G = G + (drift<1)*ml.ones((m,m))/m+(drift>1)*tau*v*ml.ones((1,m))
        
    if butools.verbose:
        if method=="PWCR":
            D = D.T
        else:
            D = Dold
        temp = D[:,-m:]
        for i in range (D.shape[1]//m-1,0,-1):
            temp = D[:,(i-1)*m:i*m] + temp*G
        res_norm = la.norm(G-temp,np.inf)
        print("Final Residual Error for G: ", res_norm)

    return G
    
def MG1StationaryDistr (A, B=None, G=None, K=500, prec=1e-14):
    """
    Returns the stationary distribution of the M/G/1 type
    Markov chain up to a given level K.
    
    Parameters
    ----------
    A : length(M) list of matrices of shape (N,N)
        Matrix blocks of the M/G/1 type generator in the 
        regular part, from 0 to M-1.
    B : length(M) list of matrices of shape (N,N)
        Matrix blocks of the M/G/1 type generator at the
        boundary, from 0 to M-1.
    G : matrix, shape (N,N)
        Matrix G of the M/G/1 type Markov chain
    K : integer
        The stationary distribution is returned up to
        this level.
    
    Returns
    -------
    pi : array, shape (1,(K+1)*N)
        The stationary probability vector up to level K
    """
    # Compute g
    if G is None:
        G = MG1FundamentalMatrix (A, prec)   
    g = DTMCSolve(G)

    A = np.hstack(A)
    m = A.shape[0]
    I = ml.eye(m)
    dega = A.shape[1]//m-1
    
    if B is None:
        mb = m
        degb = dega
        B=ml.matrix(A)
    else:
        B = np.hstack(B)
        mb = B.shape[0]
        if mb != m:
            raise Exception ("Matrix B has an incorrect number of columns")
        degb = (B.shape[1]-mb)//m
    
    # Compute theta and beta, sum_v>=0 Av and sum_v>=k Av G^v-1 
    # the last sums (for k=1,...,amax) are stored in A
    sumA = A[:,dega*m:]
    beta = np.sum(sumA, 1)
    # beta = (A_maxd)e + (A_maxd + A_maxd-1)e + ... + (Amaxd+...+A1)e
    for i in range(dega-1,0,-1):
        sumA = sumA + A[:,i*m:(i+1)*m]
        A[:,i*m:(i+1)*m] = A[:,i*m:(i+1)*m] + A[:,(i+1)*m:(i+2)*m] * G
        beta = beta + np.sum(sumA,1)

    sumA = sumA+A[:,0:m]
    theta = DTMCSolve(sumA)
    drift = (theta * beta)[0,0]
    
    if drift >= 1:
        raise Exception("The Markov chain characterized by A is not positive recurrent")
       
    if B is None:
        # Compute pi_0
        pi0 = (1.0-drift) * g
    else:
        # Compute sum_v>=1 Bv, sum_v>=1 (v-1) Bv e, sum_v>=k Bv G^v-1
        # the last sums (for k=1,...,bmax) are stored in B
        sumBB0 = ml.matrix(B[:,mb+(degb-1)*m:])
        Bbeta = np.zeros((mb,1))
        for i in range(degb-1,0,-1):
            Bbeta = Bbeta + np.sum(sumBB0,1)
            sumBB0 = sumBB0 + B[:,mb+(i-1)*m:mb+i*m]
            B[:,mb+(i-1)*m:mb+i*m] = B[:,mb+(i-1)*m:mb+i*m] + B[:,mb+i*m:mb+(i+1)*m] * G
        
        # Compute K, kappa
        Km = B[:,0:mb] + B[:,mb:mb+m] * G
        kappa = DTMCSolve(Km)

        # Compute psi1, psi2
        temp = np.sum(la.inv(I-sumA-(ml.ones((m,1))-beta)*g),1)
        psi1 = (I-A[:,0:m] - A[:,m:2*m]) * temp + 1.0/(1.0-drift) * np.sum(A[:,0:m],1)
        psi2 = ml.ones((mb,1)) + (sumBB0 - B[:,mb:mb+m]) * temp + 1.0/(1.0-drift) * Bbeta

        # Compute kappa1
        tildekappa1 = psi2 + B[:,mb:mb+m] * la.inv(I-A[:,m:2*m]) * psi1

        # Compute pi_0
        pi0 = 1.0 / (kappa*tildekappa1) * kappa
    
    numit = 1
    sumpi = np.sum(pi0)
    pi = []
    
    # Start stable RAMASWAMI formula
    invbarA1 = la.inv(I-A[:,m:2*m])
    while sumpi < 1.0-1e-10 and numit < K:
        if numit <= degb:
            if B is None:
                pix = pi0 * A[:,numit*mb:(numit+1)*mb]
            else:
                pix = pi0 * B[:,mb+(numit-1)*m:mb+numit*m]
        else:
            pix = ml.zeros((1,m))

        for j in range(1,min(numit,dega)):
            pix += pi[numit-1-j] * A[:,(j+1)*m:(j+2)*m]

        pix *= invbarA1
        sumpi = sumpi + np.sum(pix)
        pi.append (pix)
        if butools.verbose:
                print("Accumulated mass of the first ",numit, " (reblocked) components:", sumpi)
        numit+=1

    pi.insert(0, pi0)
    if butools.verbose and numit == K:
        print("Maximum Number of Components ", numit, " reached")
    return np.hstack(pi)


def GM1FundamentalMatrix (A, precision=1e-14, maxNumIt=50, method="ShiftPWCR", dual="R", maxNumRoot=2048, shiftType="one"):
    """
    Returns matrix R corresponding to the G/M/1 type Markov
    chain given by matrices A.
    
    Matrix R is the minimal non-negative solution of the 
    following matrix equation:
    
    .. math::
        R = A_0 + R A_1 + R^2 A_2 + R^3 A_3 + \dots.
    
    The implementation is based on [1]_, please cite it if
    you use this method.
    
    Parameters
    ----------
    A : length(M) list of matrices of shape (N,N)
        Matrix blocks of the G/M/1 type generator in the 
        regular part, from 0 to M-1.
    precision : double, optional
        Matrix R is computed iteratively up to this
        precision. The default value is 1e-14
    maxNumIt : int, optional
        The maximal number of iterations. The default value
        is 50.
    method : {"CR", "RR", "NI", "FI", "IS"}, optional
        The method used to solve the matrix-quadratic
        equation (CR: cyclic reduction, RR: Ramaswami
        reduction, NI: Newton iteration, FI: functional
        iteration, IS: invariant subspace method). The 
        default is "CR".
    
    Returns
    -------
    R : matrix, shape (N,N)
        The R matrix of the G/M/1 type Markov chain.
    
    References
    ----------
    .. [1] Bini, D. A., Meini, B., Steffé, S., Van Houdt,
           B. (2006, October). Structured Markov chains 
           solver: software tools. In Proceeding from the
           2006 workshop on Tools for solving structured 
           Markov chains (p. 14). ACM.
    """

    A = np.hstack(A)
    m = A.shape[0]
    I = ml.eye(m)
    dega = A.shape[1]//m-1
    
    # compute invariant vector of A and the drift
    # drift > 1: positive recurrent GIM1, drift < 1: transient GIM1
    sumA = A[:,dega*m:]
    beta = np.sum(sumA,1)
    # beta = (A_maxd)e + (A_maxd + A_maxd-1)e + ... + (Amaxd+...+A1)e
    for i in range(dega-1,0,-1):
        sumA = sumA + A[:,i*m:(i+1)*m]
        beta = beta + np.sum(sumA,1)

    sumA = sumA + A[:,:m]
    theta = DTMCSolve(sumA)
    drift = theta * beta
    
    if dual=="R" or (dual=="A" and drift <=1 ): # RAM dual
        # compute the RAM Dual process
        for i in range(dega+1):
            A[:,i*m:(i+1)*m] = Diag(1.0/theta) * A[:,i*m:(i+1)*m].T * Diag(theta)
    else: # Bright dual
        if drift > 1: # A -> positive recurrent GIM1
            # compute the Caudal characteristic of A
            eta, v = GM1TypeCaudal(A)
        else: # A -> transient GIM1 (=recurrent MG1)
            eta, v = MG1TypeDecay(A)
        # compute invariant vector of A0+A1*eta+A2*eta^2+...+Amax*eta^max
        sumAeta = eta**dega * A[:,dega*m:]
        for i in range(dega-1,-1,-1):
            sumAeta = sumAeta + eta**i * A[:,i*m:(i+1)*m]
        theta = DRPSolve(sumAeta+(1.0-eta)*I)
        # compute the Bright Dual process
        for i in range(dega+1):
            A[:,i*m:(i+1)*m] = eta**(i-1) * Diag(1.0/theta) * A[:,i*m:(i+1)*m].T * Diag(theta)
    
    G = MG1FundamentalMatrix (A, precision, maxNumIt, method, maxNumRoot, shiftType)
    
    if dual=="R" or (dual=="A" and drift <=1 ): # RAM dual
        return Diag(1.0/theta) * G.T * Diag(theta)
    else: # Bright dual
        return Diag(1.0/theta) * G.T * Diag(theta) * eta

def GM1StationaryDistr (B, R, K):
    """
    Returns the stationary distribution of the G/M/1 type
    Markov chain up to a given level K.
    
    Parameters
    ----------
    A : length(M) list of matrices of shape (N,N)
        Matrix blocks of the G/M/1 type generator in the 
        regular part, from 0 to M-1.
    B : length(M) list of matrices of shape (N,N)
        Matrix blocks of the G/M/1 type generator at the
    R : matrix, shape (N,N)
        Matrix R of the G/M/1 type Markov chain
    K : integer
        The stationary distribution is returned up to
        this level.
    
    Returns
    -------
    pi : array, shape (1,(K+1)*N)
        The stationary probability vector up to level K
    """

    if not isinstance(B,np.ndarray):
        B = np.vstack(B)

    m = R.shape[0]
    I = ml.eye(m)   

    temp = (I-R).I
    if np.max(temp<-100*butools.checkPrecision):
        raise Exception("The spectral radius of R is not below 1: GM1 is not pos. recurrent")
    
    maxb = B.shape[0]//m
    BR = B[(maxb-1)*m:,:]
    for i in range(maxb-1,0,-1):
        BR = R * BR + B[(i-1)*m:i*m,:]

    pix = DTMCSolve(BR)
    pix = pix / np.sum(pix*temp)
    
    pi = [pix]    
    sumpi = np.sum(pix)
    numit = 1
    while sumpi < 1.0-1e-10 and numit < 1+K:
        pix = pix*R; # compute pi_(numit+1)
        numit += 1
        sumpi += np.sum(pix);
        pi.append(ml.matrix(pix))
        if butools.verbose:
            print("Accumulated mass after ", numit, " iterations: ", sumpi)

    if butools.verbose and numit == K+1:
        print("Maximum Number of Components ", numit-1, " reached")
    
    return np.hstack(pi)
