import numpy as np
import numpy.matlib as ml
import scipy.linalg as la
import butools
import math
from butools.utils import Diag
from butools.reptrans import SimilarityMatrixForVectors
from butools.mam import QBDSolve
from butools.moments import MomsFromFactorialMoms
from butools.map import CheckMAPRepresentation
from butools.mc import CTMCSolve, DTMCSolve, CheckGenerator

def QBDQueue(B, L, F, L0, *argv):
    """
    Returns various performane measures of a continuous time
    QBD queue.
    
    QBD queues have a background continuous time Markov chain
    with generator Q whose the transitions can be partitioned
    into three sets: transitions accompanied by an arrival
    of a new job (F, forward), transitions accompanied by 
    the service of the current job in the server (B, 
    backward) and internal transitions (L, local). 
    Thus we have Q=B+L+F. L0 is the matrix of local 
    transition rates if the queue is empty.
    
    Parameters
    ----------
    B : matrix, shape(N,N)
        Transitions of the background process accompanied by 
        the service of the current job in the server
    L : matrix, shape(N,N)
        Internal transitions of the background process 
        that do not generate neither arrival nor service
    F : matrix, shape(N,N)
        Transitions of the background process accompanied by 
        an arrival of a new job
    L0 : matrix, shape(N,N)
        Internal transitions of the background process when
        there are no jobs in the queue
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
        | "ncDistrMG"    | None               | The vector-matrix parameters of the    |
        |                |                    | matrix-geometric distribution of the   |
        |                |                    | number of customers in the system      |
        +----------------+--------------------+----------------------------------------+
        | "ncDistrDPH"   | None               | The vector-matrix parameters of the    |
        |                |                    | matrix-geometric distribution of the   |
        |                |                    | number of customers in the system,     |
        |                |                    | converted to a discrete PH             |
        |                |                    | representation                         |
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
        |                |                    | condition when solving the             |
        |                |                    | matrix-quadratic equation              |
        +----------------+--------------------+----------------------------------------+
        
        (The quantities related to the number of customers in 
        the system include the customer in the server, and the 
        sojourn time related quantities include the service 
        times as well)
        
    Returns
    -------
    Ret : list of the performance measures
        Each entry of the list corresponds to a performance 
        measure requested. If there is just a single item, 
        then it is not put into a list.
    
    Notes
    -----
    "ncDistrMG" and "stDistrMG" behave much better numerically than 
    "ncDistrDPH" and "stDistrPH".
    """

    # parse options
    prec = 1e-14
    needST = False
    eaten = []
    for i in range(len(argv)):
        if argv[i]=="prec":
            prec = argv[i+1]
            eaten.append(i)
            eaten.append(i+1) 
        elif type(argv[i]) is str and len(argv[i])>2 and argv[i][0:2]=="st":
            needST = True
    
    if butools.checkInput and not CheckGenerator(B+L+F):
        raise Exception('QBDQueue: The matrix sum (B+L+F) is not a valid generator of a Markov chain!')
    
    if butools.checkInput and not CheckGenerator(L0+F):
        raise Exception('QBDQueue: The matrix sum (L0+F) is not a valid generator of a Markov chain!')

    pi0, R = QBDSolve (B, L, F, L0, prec)
    N = pi0.shape[1]
    I = ml.eye(N)
    
    if needST:
        U = L + R*B
        Rh = (-U).I*F
        eta = pi0*F*(I-Rh).I
        eta = eta/np.sum(eta)
        z = np.reshape(I,(N*N,1), 'F')
    
    Ret = []
    argIx = 0
    while argIx<len(argv):
        if argIx in eaten:
            argIx += 1
            continue
        elif type(argv[argIx]) is str and argv[argIx]=="ncDistrDPH":
            # transform it to DPH
            alpha = pi0*R*(I-R).I
            A = Diag(alpha).I*R.T*Diag(alpha)
            Ret.append(alpha)
            Ret.append(A)
        elif type(argv[argIx]) is str and argv[argIx]=="ncDistrMG":
            # transform it to MG
            B = SimilarityMatrixForVectors(np.sum((I-R).I*R,1), np.ones((N,1)))
            Bi = B.I
            A = B*R*Bi
            alpha = pi0*Bi
            Ret.append(alpha)
            Ret.append(A)
        elif type(argv[argIx]) is str and argv[argIx]=="ncMoms":
            numOfMoms = argv[argIx+1]
            argIx += 1
            moms = []
            iR = (I-R).I
            for m in range(1,numOfMoms+1):
                moms.append(math.factorial(m)*np.sum(pi0*iR**(m+1)*R**m))
            Ret.append(MomsFromFactorialMoms(moms))
        elif type(argv[argIx]) is str and argv[argIx]=="ncDistr":
            numOfQLProbs = argv[argIx+1]
            argIx += 1
            values = np.empty(numOfQLProbs)
            values[0] = np.sum(pi0)
            RPow = I
            for p in range(numOfQLProbs-1):
                RPow = RPow * R
                values[p+1] = np.sum(pi0*RPow)
            Ret.append(values)
        elif type(argv[argIx]) is str and argv[argIx]=="stDistrPH":
            # transform to ph distribution
            ix = np.arange(N)
            nz = ix[eta.flat>prec]
            Delta = Diag(eta)
            A = np.kron(L+F,I[nz,:][:,nz]) + np.kron(B,Delta[nz,:][:,nz].I*Rh[nz,:][:,nz].T*Delta[nz,:][:,nz])
            alpha = z.T*np.kron(I,Delta[:,nz])
            Ret.append(alpha)
            Ret.append(A)
        elif type(argv[argIx]) is str and argv[argIx]=="stDistrME":
            # transform it such that the closing vector is a vector of ones
            # this is the way butools accepts ME distributions
            Bm = SimilarityMatrixForVectors(z, np.ones(z.shape))
            Bmi = Bm.I
            A = Bm * (np.kron(L.T+F.T,I) + np.kron(B.T,Rh)) * Bmi
            alpha = np.kron(ml.ones((1,N)), eta) * Bmi
            Ret.append(alpha)
            Ret.append(A)
        elif type(argv[argIx]) is str and argv[argIx]=="stMoms":
            numOfMoms = argv[argIx+1]
            argIx += 1;
            moms = []
            Z = np.kron(L.T+F.T,I)+np.kron(B.T,Rh)
            iZ = -Z.I
            for m in range(1,numOfMoms+1):
                moms.append(math.factorial(m)*np.sum(np.kron(ml.ones((1,N)), eta)*iZ**(m+1)*(-Z)*z))
            Ret.append(moms)
        elif type(argv[argIx]) is str and argv[argIx]=="stDistr":
            points = argv[argIx+1]
            argIx += 1
            values = np.empty(points.shape)
            Z = np.kron(L.T+F.T,I)+np.kron(B.T,Rh)
            for p in range(len(points.flat)):
                values.flat[p] = 1.0-np.sum(np.kron(ml.ones((1,N)), eta)*la.expm(Z*points.flat[p])*z)
            Ret.append(values)
        else:
            raise Exception("QBDQueue: Unknown parameter "+str(argv[argIx]))
        argIx += 1
    if len(Ret)==1:
        return Ret[0]
    else:
        return Ret

def MAPMAP1(D0, D1, S0, S1, *argv):
    """
    Returns various performane measures of a continuous time
    MAP/MAP/1 queue.
    
    In a MAP/MAP/1 queue both the arrival and the service
    processes are characterized by Markovian arrival 
    processes.
    
    Parameters
    ----------
    D0 : matrix, shape(N,N)
        The transitions of the arrival MAP not accompanied by
        job arrivals
    D1 : matrix, shape(N,N)
        The transitions of the arrival MAP accompanied by
        job arrivals
    S0 : matrix, shape(N,N)
        The transitions of the service MAP not accompanied by
        job service
    S1 : matrix, shape(N,N)
        The transitions of the service MAP accompanied by
        job service
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
        | "ncDistrMG"    | None               | The vector-matrix parameters of the    |
        |                |                    | matrix-geometric distribution of the   |
        |                |                    | number of customers in the system      |
        +----------------+--------------------+----------------------------------------+
        | "ncDistrDPH"   | None               | The vector-matrix parameters of the    |
        |                |                    | matrix-geometric distribution of the   |
        |                |                    | number of customers in the system,     |
        |                |                    | converted to a discrete PH             |
        |                |                    | representation                         |
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
        |                |                    | condition when solving the             |
        |                |                    | matrix-quadratic equation              |
        +----------------+--------------------+----------------------------------------+
        
        (The quantities related to the number of customers in 
        the system include the customer in the server, and the 
        sojourn time related quantities include the service 
        times as well)
    
    Returns
    -------
    Ret : list of the performance measures
        Each entry of the list corresponds to a performance 
        measure requested. If there is just a single item, 
        then it is not put into a list.
        Notes
    -----
    "ncDistrMG" and "stDistrME" behave much better numerically than 
    "ncDistrDPH" and "stDistrPH".
    """

    # parse options
    prec = 1e-14
    needST = False
    eaten = []
    for i in range(len(argv)):
        if argv[i]=="prec":
            prec = argv[i+1]
            eaten.append(i)
            eaten.append(i+1) 
        elif type(argv[i]) is str and len(argv[i])>2 and argv[i][0:2]=="st":
            needST = True
    
    if butools.checkInput and not CheckMAPRepresentation(D0,D1):
        raise Exception('MAPMAP1: The arrival process (D0,D1) is not a valid MAP representation!')
    
    if butools.checkInput and not CheckMAPRepresentation(S0,S1):
        raise Exception('MAPMAP1: The service process (S0,S1) is not a valid MAP representation!')

    IA = ml.eye(D0.shape[0])
    IS = ml.eye(S0.shape[0])
    
    B = np.kron(IA,S1)
    L = np.kron(D0,IS)+np.kron(IA,S0)
    F = np.kron(D1,IS)
    L0 = np.kron(D0,IS)
    
    pi0, R = QBDSolve (B, L, F, L0, prec)
    N = pi0.shape[1]
    I = ml.eye(N)
    
    if needST:
        # calculate the distribution of the age at departures
        U = L + R*B
        Rh = -U.I*F
        T = np.kron(IA,S0) + Rh * B
        eta = pi0*F*(I-Rh).I
        eta = eta/np.sum(eta)
    
    Ret = []
    argIx = 0
    while argIx<len(argv):
        if argIx in eaten:
            argIx += 1
            continue
        elif type(argv[argIx]) is str and argv[argIx]=="ncDistrDPH":
            # transform it to DPH
            alpha = pi0*R*(I-R).I
            A = Diag(alpha).I*R.T*Diag(alpha)
            Ret.append(alpha)
            Ret.append(A)
        elif type(argv[argIx]) is str and argv[argIx]=="ncDistrMG":
            # transform it to MG
            B = SimilarityMatrixForVectors(np.sum((I-R).I*R,1), np.ones((N,1)))
            Bi = B.I
            A = B*R*Bi
            alpha = pi0*Bi
            Ret.append(alpha)
            Ret.append(A)
        elif type(argv[argIx]) is str and argv[argIx]=="ncMoms":
            numOfMoms = argv[argIx+1]
            argIx += 1
            moms = []
            iR = (I-R).I
            for m in range(1,numOfMoms+1):
                moms.append(math.factorial(m)*np.sum(pi0*iR**(m+1)*R**m))
            Ret.append(MomsFromFactorialMoms(moms))
        elif type(argv[argIx]) is str and argv[argIx]=="ncDistr":
            numOfQLProbs = argv[argIx+1]
            argIx += 1
            values = np.empty(numOfQLProbs)
            values[0] = np.sum(pi0)
            RPow = I
            for p in range(numOfQLProbs-1):
                RPow = RPow * R
                values[p+1] = np.sum(pi0*RPow)
            Ret.append(values)
        elif type(argv[argIx]) is str and argv[argIx]=="stDistrPH":
            # transform it to PH representation
            beta = CTMCSolve(S0+S1)
            theta = DTMCSolve(-D0.I*D1)
            vv = np.kron(theta,beta)
            ix = np.arange(N)
            nz = ix[vv.flat>prec]
            delta = Diag(vv[:,nz])
            alpha = ml.ones((1,N))*B[nz,:].T*delta / np.sum(beta*S1)
            A = delta.I*T[nz,:][:,nz].T*delta
            Ret.append(alpha)
            Ret.append(A)
        elif type(argv[argIx]) is str and argv[argIx]=="stDistrME":
            Ret.append(eta)
            Ret.append(T)
        elif type(argv[argIx]) is str and argv[argIx]=="stMoms":
            numOfMoms = argv[argIx+1]
            argIx += 1
            moms = []
            iT = -T.I
            for m in range(1,numOfMoms+1):
                moms.append(math.factorial(m)*np.sum(eta*iT**m))
            Ret.append(moms)
        elif type(argv[argIx]) is str and argv[argIx]=="stDistr":
            points = argv[argIx+1]
            argIx += 1
            values = np.empty(points.shape)
            for p in range(len(points.flat)):
                values.flat[p] = 1.0-np.sum(eta*la.expm(T*points.flat[p]))
            Ret.append(values)
        else:
            raise Exception("MAPMAP1: Unknown parameter "+str(argv[argIx]))
        argIx += 1
    if len(Ret)==1:
        return Ret[0]
    else:
        return Ret
