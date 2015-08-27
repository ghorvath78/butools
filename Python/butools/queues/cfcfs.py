import numpy as np
import numpy.matlib as ml
import scipy.linalg as la
import butools
import math
from butools.mam import GeneralFluidSolve
from butools.utils import Linsolve, Diag
from butools.reptrans import SimilarityMatrixForVectors
from butools.mc import CTMCSolve, CheckGenerator

def FluidQueue (Q, Rin, Rout, *argv):
    """
    Returns various performane measures of a fluid queue.
    
    In a fluid queue there is a background continuous time
    Markov chain (given by generator Q), and diagonal
    matrix Rin (Rout) whose ith entry provides the 
    fluid rate at which fluid enters the queue (can be 
    served) while the background process is in state i.
    
    Parameters
    ----------
    Q : matrix, shape (N,N)
        The generator of the background Markov chain
    Rin : matrix, shape (N,N)
        Diagonal matrix containing the fluid input rates
        associated to the states of the background process
    Rout : matrix, shape (N,N)
        Diagonal matrix containing the fluid output rates
        associated to the states of the background process
    further parameters : 
        The rest of the function parameters specify the options
        and the performance measures to be computed.
    
        The supported performance measures and options in this 
        function are:
    
        +----------------+--------------------+--------------------------------------+
        | Parameter name | Input parameters   | Output                               |
        +================+====================+======================================+
        | "flMoms"       | Number of moments  | The moments of the fluid level       |
        +----------------+--------------------+--------------------------------------+
        | "flDistr"      | A vector of points | The fluid level distribution at      |
        |                |                    | the requested points (cdf)           |
        +----------------+--------------------+--------------------------------------+
        | "flDistrME"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-exponentially distributed     |
        |                |                    | fluid level distribution             |
        +----------------+--------------------+--------------------------------------+
        | "flDistrPH"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-exponentially distributed     |
        |                |                    | fluid level distribution, converted  |
        |                |                    | to a PH representation               |
        +----------------+--------------------+--------------------------------------+
        | "stMoms"       | Number of moments  | The sojourn time moments of fluid    |
        |                |                    | drops                                |
        +----------------+--------------------+--------------------------------------+
        | "stDistr"      | A vector of points | The sojourn time distribution at the |
        |                |                    | requested points (cummulative, cdf)  |
        +----------------+--------------------+--------------------------------------+
        | "stDistrME"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-exponentially distributed     |
        |                |                    | sojourn time distribution            |
        +----------------+--------------------+--------------------------------------+
        | "stDistrPH"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-exponentially distributed     |
        |                |                    | sojourn time distribution, converted |
        |                |                    | to a PH representation               |
        +----------------+--------------------+--------------------------------------+
        | "prec"         | The precision      | Numerical precision to check if the  |
        |                |                    | input is valid and it is also used   |
        |                |                    | as a stopping condition when solving |
        |                |                    | the Riccati equation                 |
        +----------------+--------------------+--------------------------------------+
        | "Q0"           | Matrix, shape(N,N) | The generator of the background      |
        |                |                    | Markov chain when the fluid level is |
        |                |                    | zero. If not given, Q0=Q is assumed  |
        +----------------+--------------------+--------------------------------------+
    
    Returns
    -------
    Ret : list of the performance measures
        Each entry of the list corresponds to a performance 
        measure requested. If there is just a single item, 
        then it is not put into a list.
    
    Notes
    -----
    "flDistrME" and "stDistrME" behave much better numerically than 
    "flDistrPH" and "stDistrPH".
    """

    # parse options
    prec = 1e-14
    needST = False
    Q0 = []
    eaten = []
    for i in range(len(argv)):
        if type(argv[i]) is str and argv[i]=="prec":
            prec = argv[i+1]
            eaten.append(i)
            eaten.append(i+1)
        elif type(argv[i]) is str and argv[i]=="Q0":
            Q0 = argv[i+1]
            eaten.append(i)
            eaten.append(i+1)
        elif type(argv[i]) is str and len(argv[i])>2 and argv[i][0:2]=="st":
            needST = True

    if butools.checkInput and not CheckGenerator(Q,False):
        raise Exception('FluidQueue: Generator matrix Q is not Markovian!')

    if butools.checkInput and len(Q0)>0 and not CheckGenerator(Q0,False):
        raise Exception('FluidQueue: Generator matrix Q0 is not Markovian!')

    if butools.checkInput and (np.any(np.diag(Rin)<-butools.checkPrecision) or np.any(np.diag(Rout)<-butools.checkPrecision)):
        raise Exception('FluidQueue: Fluid rates Rin and Rout must be non-negative !')

    mass0, ini, K, clo = GeneralFluidSolve (Q, Rin-Rout, Q0, prec)
    if needST:
        N = Q.shape[0]
        iniKi = Linsolve(K.T,-ini.T).T # iniki = ini*inv(-K);
        lambd = np.sum(mass0*Rin + iniKi*clo*Rin)
    
    Ret = []
    argIx = 0
    while argIx<len(argv):
        if argIx in eaten:
            argIx += 1
            continue
        elif type(argv[argIx]) is str and argv[argIx]=='flDistrPH':
            # transform it to PH
            Delta = Diag(Linsolve(K.T,-ini.T)) # Delta = diag (ini*inv(-K));
            A = Delta.I*K.T*Delta
            alpha = np.sum(clo,1).T*Delta
            Ret.append(alpha)
            Ret.append(A)
        elif type(argv[argIx]) is str and argv[argIx]=='flDistrME':
            # transform it to ME
            B = SimilarityMatrixForVectors((-K).I*np.sum(clo,1), np.ones((K.shape[0],1)))
            Bi = B.I
            alpha = ini*Bi
            A = B*K*Bi
            Ret.append(alpha)
            Ret.append(A)
        elif type(argv[argIx]) is str and argv[argIx]=='flMoms':
            numOfMoms = argv[argIx+1]
            argIx += 1
            moms = []
            iK = -K.I
            for m in range(1,numOfMoms+1):
                moms.append(math.factorial(m)*np.sum(ini*iK**(m+1)*clo))
            Ret.append(moms)
        elif type(argv[argIx]) is str and argv[argIx]=='flDistr':
            points = argv[argIx+1]
            argIx += 1
            values = np.empty(points.shape)
            iK = -K.I
            for p in range(len(points.flat)):
                values.flat[p] = np.sum(mass0) + np.sum(ini*(ml.eye(K.shape[0])-la.expm(K*points.flat[p]))*iK*clo)
            Ret.append (values)
        elif type(argv[argIx]) is str and argv[argIx]=='stDistrPH':
            # transform it to PH
            Delta = Diag(iniKi/lambd)
            alpha = np.reshape(clo*Rin,(1,N*len(ini.flat)),'F')*np.kron(ml.eye(N),Delta)
            A = np.kron(Rout, Delta.I*K.T*Delta) + np.kron(Q, ml.eye(K.shape[0]))
            Ret.append(alpha)
            Ret.append(A)
        elif type(argv[argIx]) is str and argv[argIx]=='stDistrME':
            B = SimilarityMatrixForVectors(np.reshape(-K.I*clo*Rin,(N*ini.size,1),'F'), np.ones((N*ini.size,1)))
            Bi = B.I
            alpha = np.kron(ml.ones((1,N)), ini/lambd)*Bi
            A = B*(np.kron(Q.T,ml.eye(K.shape[0])) + np.kron(Rout,K))*Bi        
            Ret.append(alpha)
            Ret.append(A)
        elif type(argv[argIx]) is str and argv[argIx]=='stMoms':
            numOfMoms = argv[argIx+1]
            argIx += 1
            moms = []
            Z = np.kron(Q.T,ml.eye(K.shape[0])) + np.kron(Rout,K)
            iZ = -Z.I
            kini = np.kron(ml.ones((1,N)), ini/lambd)
            kclo = np.reshape(-K.I*clo*Rin,(N*ini.size,1),'F')
            for m in range(1,numOfMoms+1):
                moms.append(math.factorial(m)*np.sum(kini*iZ**(m+1)*(-Z)*kclo))
            Ret.append(moms)
        elif type(argv[argIx]) is str and argv[argIx]=='stDistr':
            points = argv[argIx+1]
            argIx += 1
            values = np.empty(points.shape)
            Z = np.kron(Q.T,ml.eye(K.shape[0])) + np.kron(Rout,K)
            kini = np.kron(ml.ones((1,N)), ini/lambd)
            kclo = np.reshape(-K.I*clo*Rin,(N*ini.size,1),'F')
            for p in range(len(points.flat)):
                values.flat[p] = 1.0-np.sum(kini*la.expm(Z*points.flat[p])*kclo)
            Ret.append(values)
        else:
            raise Exception ("FluidQueue: Unknown parameter "+str(argv[argIx]))
        argIx += 1

    if len(Ret)==1:
        return Ret[0]
    else:
        return Ret

def FluFluQueue(Qin, Rin, Qout, Rout, srv0stop, *argv):
    """
    Returns various performane measures of a fluid queue
    with independent fluid arrival and service processes.
    
    Two types of boundary behavior is available. If 
    srv0stop=false, the output process evolves continuously
    even if the queue is empty. If srv0stop=true, the 
    output process slows down if there is fewer fluid in
    the queue than it can serve. If the queue is empty
    and the fluid input rate is zero, the output process
    freezes till fluid arrives.
    
    Parameters
    ----------
    Qin : matrix, shape (N,N)
        The generator of the background Markov chain 
        corresponding to the input process
    Rin : matrix, shape (N,N)
        Diagonal matrix containing the fluid input rates
        associated to the states of the input background 
        process
    Qout : matrix, shape (N,N)
        The generator of the background Markov chain 
        corresponding to the output process
    Rout : matrix, shape (N,N)
        Diagonal matrix containing the fluid output rates
        associated to the states of the input background 
        process
    srv0stop : bool
        If true, the service output process slows down if
        there is fewer fluid in the queue than it can 
        serve. If false, the output process evolves 
        continuously.
    further parameters : 
        The rest of the function parameters specify the options
        and the performance measures to be computed.
    
        The supported performance measures and options in this 
        function are:
        
        +----------------+--------------------+--------------------------------------+
        | Parameter name | Input parameters   | Output                               |
        +================+====================+======================================+
        | "flMoms"       | Number of moments  | The moments of the fluid level       |
        +----------------+--------------------+--------------------------------------+
        | "flDistr"      | A vector of points | The fluid level distribution at      |
        |                |                    | the requested points (cdf)           |
        +----------------+--------------------+--------------------------------------+
        | "flDistrME"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-exponentially distributed     |
        |                |                    | fluid level distribution             |
        +----------------+--------------------+--------------------------------------+
        | "flDistrPH"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-exponentially distributed     |
        |                |                    | fluid level distribution, converted  |
        |                |                    | to a PH representation               |
        +----------------+--------------------+--------------------------------------+
        | "stMoms"       | Number of moments  | The sojourn time moments of fluid    |
        |                |                    | drops                                |
        +----------------+--------------------+--------------------------------------+
        | "stDistr"      | A vector of points | The sojourn time distribution at the |
        |                |                    | requested points (cummulative, cdf)  |
        +----------------+--------------------+--------------------------------------+
        | "stDistrME"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-exponentially distributed     |
        |                |                    | sojourn time distribution            |
        +----------------+--------------------+--------------------------------------+
        | "stDistrPH"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-exponentially distributed     |
        |                |                    | sojourn time distribution, converted |
        |                |                    | to a PH representation               |
        +----------------+--------------------+--------------------------------------+
        | "prec"         | The precision      | Numerical precision to check if the  |
        |                |                    | input is valid and it is also used   |
        |                |                    | as a stopping condition when solving |
        |                |                    | the Riccati equation                 |
        +----------------+--------------------+--------------------------------------+
    
    Returns
    -------
    Ret : list of the performance measures
        Each entry of the list corresponds to a performance 
        measure requested. If there is just a single item, 
        then it is not put into a list.
    
    Notes
    -----
    "flDistrME" and "stDistrME" behave much better numerically than 
    "flDistrPH" and "stDistrPH".
    
    References
    ----------
    .. [1] Horvath G, Telek M, "Sojourn times in fluid queues 
           with independent and dependent input and output 
           processes PERFORMANCE EVALUATION 79: pp. 160-181, 2014.
    """

    # parse options
    prec = 1e-14
    needST = False
    needQL = False
    Q0 = []
    eaten = []
    for i in range(len(argv)):
        if type(argv[i]) is str and argv[i]=="prec":
            prec = argv[i+1]
            eaten.append(i)
            eaten.append(i+1)
        elif type(argv[i]) is str and len(argv[i])>2 and argv[i][0:2]=="st":
            needST = True
        elif type(argv[i]) is str and len(argv[i])>2 and argv[i][0:2]=="fl":
            needQL = True

    if butools.checkInput and not CheckGenerator(Qin,False):
        raise Exception('FluFluQueue: Generator matrix Qin is not Markovian!')

    if butools.checkInput and not CheckGenerator(Qout,False):
        raise Exception('FluFluQueue: Generator matrix Qout is not Markovian!')

    if butools.checkInput and (np.any(np.diag(Rin)<-butools.checkPrecision) or np.any(np.diag(Rout)<-butools.checkPrecision)):
        raise Exception('FluFluQueue: Fluid rates Rin and Rout must be non-negative !')

    Iin = ml.eye(Qin.shape[0])
    Iout = ml.eye(Qout.shape[0])

    if needQL:
        Q = np.kron(Qin,Iout)+np.kron(Iin,Qout)
        if srv0stop:
            Q0 = np.kron(Qin,Iout)+np.kron(Rin, la.pinv(Rout)*Qout)
        else:
            Q0 = Q
        mass0, ini, K, clo = GeneralFluidSolve (Q, np.kron(Rin,Iout)-np.kron(Iin,Rout), Q0, prec)

    if needST:
        Rh = np.kron(Rin,Iout) - np.kron(Iin,Rout)
        Qh = np.kron(Qin, Rout) + np.kron(Rin, Qout)
        massh, inih, Kh, cloh = GeneralFluidSolve (Qh, Rh, prec=prec)

        # sojourn time density in case of 
        # srv0stop = false: inih*expm(Kh*x)*cloh*kron(Rin,Iout)/lambda
        # srv0stop = true: inih*expm(Kh*x)*cloh*kron(Rin,Rout)/lambda/mu    
        lambd = np.sum(CTMCSolve(Qin)*Rin)
        mu = np.sum(CTMCSolve(Qout)*Rout)
    
    Ret = []
    argIx = 0
    while argIx<len(argv):
        if argIx in eaten:
            argIx += 1
            continue
        elif type(argv[argIx]) is str and argv[argIx]=='flDistrPH':
            # transform it to PH
            Delta = Diag(Linsolve(K.T,-ini.T)) # Delta = diag (ini*inv(-K));
            A = Delta.I*K.T*Delta
            alpha = np.sum(clo,1).T*Delta
            Ret.append(alpha)
            Ret.append(A)
        elif type(argv[argIx]) is str and argv[argIx]=='flDistrME':
            # transform it to ME
            B = SimilarityMatrixForVectors((-K).I*np.sum(clo,1), np.ones((K.shape[0],1)))
            Bi = B.I
            alpha = ini*Bi
            A = B*K*Bi
            Ret.append(alpha)
            Ret.append(A)
        elif type(argv[argIx]) is str and argv[argIx]=='flMoms':
            numOfMoms = argv[argIx+1]
            argIx += 1
            moms = []
            iK = -K.I
            for m in range(1,numOfMoms+1):
                moms.append(math.factorial(m)*np.sum(ini*iK**(m+1)*clo))
            Ret.append(moms)
        elif type(argv[argIx]) is str and argv[argIx]=='flDistr':
            points = argv[argIx+1]
            argIx += 1
            values = np.empty(points.shape)
            iK = -K.I
            for p in range(len(points.flat)):
                values.flat[p] = np.sum(mass0) + np.sum(ini*(ml.eye(K.shape[0])-la.expm(K*points.flat[p]))*iK*clo)
            Ret.append (values)
        elif type(argv[argIx]) is str and argv[argIx]=='stDistrPH':
            # convert result to PH representation
            Delta = Diag(Linsolve(Kh.T,-inih.T)) # Delta = diag (inih*inv(-Kh));
            A = Delta.I*Kh.T*Delta
            if not srv0stop:
                alpha = np.sum(Delta*cloh*np.kron(Rin,Iout)/lambd,1).T
            else:
                alpha = np.sum(Delta*cloh*np.kron(Rin,Rout)/lambd/mu,1).T
            Ret.append(alpha)
            Ret.append(A)
        elif type(argv[argIx]) is str and argv[argIx]=='stDistrME':
            # convert result to ME representation
            if not srv0stop:
                B = SimilarityMatrixForVectors(np.sum(cloh*np.kron(Rin,Iout)/lambd,1), np.ones((Kh.shape[0],1)))
            else:
                B = SimilarityMatrixForVectors(np.sum(cloh*np.kron(Rin,Rout)/lambd/mu,1), np.ones((Kh.shape[0],1)))
            iB = B.I
            A = B*Kh*iB
            alpha = inih*(-Kh).I*iB
            Ret.append(alpha)
            Ret.append(A)
        elif type(argv[argIx]) is str and argv[argIx]=='stMoms':
            numOfMoms = argv[argIx+1]
            argIx += 1
            moms = []
            if srv0stop:
                kclo = cloh*np.kron(Rin,Rout)/lambd/mu
            else:
                kclo = cloh*np.kron(Rin,Iout)/lambd
            iKh = -Kh.I
            for m in range(1,numOfMoms+1):
                moms.append(math.factorial(m)*np.sum(inih*iKh**(m+1)*kclo))
            Ret.append(moms)
        elif type(argv[argIx]) is str and argv[argIx]=='stDistr':
            points = argv[argIx+1]
            argIx += 1
            values = np.empty(points.shape)
            if srv0stop:
                kclo = cloh*np.kron(Rin,Rout)/lambd/mu
            else:
                kclo = cloh*np.kron(Rin,Iout)/lambd
            iKh = -Kh.I
            for p in range(len(points.flat)):
                values.flat[p] = 1.0-np.sum(inih*la.expm(Kh*points.flat[p])*iKh*kclo)
            Ret.append(values)
        else:
            raise Exception ("FluFluQueue: Unknown parameter "+str(argv[argIx]))
        argIx += 1

    if len(Ret)==1:
        return Ret[0]
    else:
        return Ret
