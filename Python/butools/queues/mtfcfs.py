import numpy as np
import numpy.matlib as ml
import scipy.linalg as la
import butools
import math
from butools.moments import MomsFromFactorialMoms
from butools.map import CheckMMAPRepresentation
from butools.mc import CTMCSolve
from butools.ph import CheckPHRepresentation
from butools.mam import FluidFundamentalMatrices
from butools.utils import Diag
from butools.reptrans import SimilarityMatrixForVectors

def MMAPPH1FCFS(D, sigma, S, *argv):
    """
    Returns various performane measures of a MMAP[K]/PH[K]/1 
    first-come-first-serve queue, see [1]_.
    
    Parameters
    ----------
    D : list of matrices of shape (N,N), length (K+1)
        The D0...DK matrices of the arrival process.
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
        | "stDistrME"    | None               | The vector-matrix parameters of the    |
        |                |                    | matrix-exponentially distributed       |
        |                |                    | sojourn time distribution              |
        +----------------+--------------------+----------------------------------------+
        | "stDistrPH"    | None               | The vector-matrix parameters of the    |
        |                |                    | matrix-exponentially distributed       |
        |                |                    | sojourn time distribution, converted   |
        |                |                    | to a continuous PH representation      |
        +----------------+--------------------+----------------------------------------+
        | "prec"         | The precision      | Numerical precision used as a stopping |
        |                |                    | condition when solving the Riccati     |
        |                |                    | equation                               |
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
    .. [1] Qiming He, "Analysis of a continuous time 
           SM[K]/PH[K]/1/FCFS queue: Age process, sojourn times,
           and queue lengths", Journal of Systems Science and 
           Complexity, 25(1), pp 133-155, 2012.
    """
    
    K = len(D)-1

    # parse options
    eaten = []
    precision = 1e-14;
    classes = np.arange(0,K)
    for i in range(len(argv)):
        if argv[i]=="prec":
            precision = argv[i+1]
            eaten.append(i)
            eaten.append(i+1) 
        elif argv[i]=="classes":
            classes = np.array(argv[i+1])-1
            eaten.append(i)
            eaten.append(i+1) 
    
    if butools.checkInput and not CheckMMAPRepresentation(D):
        raise Exception('MMAPPH1FCFS: The arrival process is not a valid MMAP representation!')
    
    if butools.checkInput:
        for k in range(K):
            if not CheckPHRepresentation(sigma[k],S[k]):
                raise Exception('MMAPPH1FCFS: the vector and matrix describing the service times is not a valid PH representation!')

    # some preparation
    D0 = D[0]
    N = D0.shape[0]
    Ia = ml.eye(N);
    Da = ml.zeros((N,N))
    for q in range(K):
        Da += D[q+1]
    theta = CTMCSolve(D0+Da)
    beta = [CTMCSolve(S[k]+ml.sum(-S[k],1)*sigma[k]) for k in range(K)]
    lambd = [np.sum(theta*D[k+1]) for k in range(K)]    
    mu = [np.sum(beta[k]*(-S[k])) for k in range(K)]
    Nsk = [S[k].shape[0] for k in range(K)]    
    ro = np.sum(np.array(lambd)/np.array(mu))
    alpha = theta*Da/sum(lambd)
    D0i = (-D0).I

    Sa = S[0];
    sa = [ml.zeros(sigma[0].shape)]*K
    sa[0] = sigma[0]
    ba = [ml.zeros(beta[0].shape)]*K
    ba[0] = beta[0]
    sv = [ml.zeros((Nsk[0],1))]*K
    sv[0] = ml.sum(-S[0],1)
    Pk = [D0i*D[q+1] for q in range(K)]

    for k in range(1,K):
        Sa = la.block_diag(Sa, S[k])
        for q in range(K):
            if q==k:
                sa[q] = ml.hstack((sa[q], sigma[k]))
                ba[q] = ml.hstack((ba[q], beta[k]))
                sv[q] = ml.vstack((sv[q], -np.sum(S[k],1)))
            else:
                sa[q] = ml.hstack((sa[q], ml.zeros(sigma[k].shape)))
                ba[q] = ml.hstack((ba[q], ml.zeros(beta[k].shape)))
                sv[q] = ml.vstack((sv[q], ml.zeros((Nsk[k],1))))
    Sa = ml.matrix(Sa)
    P = D0i*Da
    iVec = ml.kron(D[1],sa[0])
    for k in range(1,K):
        iVec += ml.kron(D[k+1],sa[k])
    Ns = Sa.shape[0]
    Is = ml.eye(Ns)
    
    # step 1. solve the age process of the queue
    # ==========================================

    # solve Y0 and calculate T
    Y0 = FluidFundamentalMatrices (ml.kron(Ia,Sa), ml.kron(Ia,-ml.sum(Sa,1)), iVec, D0, "P", precision)
    T = ml.kron(Ia,Sa) + Y0 * iVec
    
    # calculate pi0 and v0
    pi0 = ml.zeros((1,T.shape[0]))
    for k in range(K):
        pi0 += ml.kron(theta*D[k+1],ba[k]/mu[k])
    pi0 = - pi0 * T

    iT = (-T).I
    oa = ml.ones((N,1))

    # step 2. calculate performance measures
    # ======================================
    Ret = []
    for k in classes:
        argIx = 0
        clo = iT*ml.kron(oa,sv[k])
        while argIx<len(argv):
            if argIx in eaten:
                argIx += 1
                continue
            elif type(argv[argIx]) is str and argv[argIx]=="stMoms":
                numOfSTMoms = argv[argIx+1]
                rtMoms = []
                for m in range(1,numOfSTMoms+1):
                    rtMoms.append(math.factorial(m) * np.sum(pi0 * iT**m * clo / (pi0*clo)))
                Ret.append(rtMoms)
                argIx += 1
            elif type(argv[argIx]) is str and argv[argIx]=="stDistr":
                stCdfPoints = argv[argIx+1]
                cdf = [];
                for t in stCdfPoints:
                    pr = 1 - np.sum(pi0 * la.expm(T*t) * clo / (pi0*clo))
                    cdf.append(pr)
                Ret.append(np.array(cdf))
                argIx += 1
            elif type(argv[argIx]) is str and argv[argIx]=="stDistrME":
                Bm = SimilarityMatrixForVectors(clo/(pi0*clo),ml.ones((N*Ns,1)))
                Bmi = Bm.I
                A = Bm * T * Bmi
                alpha = pi0 * Bmi
                Ret.append(alpha)
                Ret.append(A)
            elif type(argv[argIx]) is str and argv[argIx]=="stDistrPH":
                vv = pi0*iT
                ix = np.arange(N*Ns)
                nz = ix[vv.flat>precision]
                delta = Diag(vv[:,nz])
                cl = -T*clo/(pi0*clo)
                alpha = cl[nz,:].T*delta
                A = delta.I*T[nz,:][:,nz].T*delta
                Ret.append(alpha)
                Ret.append(A)
            elif type(argv[argIx]) is str and argv[argIx]=="ncDistr":
                numOfQLProbs = argv[argIx+1]
                argIx += 1
                values = np.empty(numOfQLProbs)
                jm = ml.zeros((Ns,1))
                jm[np.sum(Nsk[0:k]):np.sum(Nsk[0:k+1]),:] = 1
                jmc = ml.ones((Ns,1))
                jmc[np.sum(Nsk[0:k]):np.sum(Nsk[0:k+1]),:] = 0
                LmCurr = la.solve_sylvester(T, ml.kron(D0+Da-D[k+1],Is), -ml.eye(N*Ns))
                values[0] = 1-ro+np.sum(pi0*LmCurr*ml.kron(oa,jmc))
                for i in range(1,numOfQLProbs):
                    LmPrev = LmCurr
                    LmCurr = la.solve_sylvester(T, ml.kron(D0+Da-D[k+1],Is), -LmPrev*ml.kron(D[k+1],Is))
                    values[i] = np.sum(pi0*LmCurr*ml.kron(oa,jmc) + pi0*LmPrev*ml.kron(oa,jm));
                Ret.append(values)
            elif type(argv[argIx]) is str and argv[argIx]=="ncMoms":
                numOfQLMoms = argv[argIx+1]
                argIx += 1
                jm = ml.zeros((Ns,1))
                jm[np.sum(Nsk[0:k]):np.sum(Nsk[0:k+1]),:] = 1
                ELn = [la.solve_sylvester(T, ml.kron(D0+Da,Is), -ml.eye(N*Ns))]
                qlMoms = []
                for n in range(1,numOfQLMoms+1):
                    bino = 1
                    Btag = ml.zeros((N*Ns,N*Ns))
                    for i in range(n):
                        Btag += bino * ELn[i]
                        bino *= (n-i) / (i+1)
                    ELn.append(la.solve_sylvester(T, ml.kron(D0+Da,Is), -Btag*ml.kron(D[k+1],Is)))
                    qlMoms.append(np.sum(pi0*ELn[n]) + np.sum(pi0*Btag*ml.kron(oa,jm)))
                Ret.append(qlMoms)
            else:
                raise Exception("MMAPPH1FCFS: Unknown parameter "+str(argv[argIx]))
            argIx += 1

    if len(Ret)==1:
        return Ret[0]
    else:
        return Ret
        
