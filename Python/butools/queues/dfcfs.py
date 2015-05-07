import numpy as np
import numpy.matlib as ml
import scipy.linalg as la
import butools
import math
from butools.utils import Diag
from butools.reptrans import TransformToOnes
from butools.mam import QBDSolve
from butools.moments import MomsFromFactorialMoms
from butools.map import CheckMAPRepresentation
from butools.mc import CTMCSolve, DTMCSolve, CheckGenerator

def QBDQueue(B, L, F, L0, *argv):

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
    
    if butools.checkInput and not CheckGenerator(B+L+F,prec):
        raise Exception('QBDQueue: The matrix sum (B+L+F) is not a valid generator of a Markov chain!')
    
    if butools.checkInput and not CheckGenerator(L0+F,prec):
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
        elif type(argv[argIx]) is str and argv[argIx]=="qlDistrDPH":
            # transform it to DPH
            alpha = pi0*R*(I-R).I
            A = Diag(alpha).I*R.T*Diag(alpha)
            Ret.append((alpha, A))
        elif type(argv[argIx]) is str and argv[argIx]=="qlDistrMG":
            # transform it to MG
            B = TransformToOnes(np.sum((I-R).I*R,1))
            Bi = B.I
            A = B*R*Bi
            alpha = pi0*Bi
            Ret.append((alpha, A))
        elif type(argv[argIx]) is str and argv[argIx]=="qlMoms":
            numOfMoms = argv[argIx+1]
            argIx += 1
            moms = []
            iR = (I-R).I
            for m in range(1,numOfMoms+1):
                moms.append(math.factorial(m)*np.sum(pi0*iR**(m+1)*R**m))
            Ret.append(MomsFromFactorialMoms(moms))
        elif type(argv[argIx]) is str and argv[argIx]=="qlDistr":
            points = argv[argIx+1]
            argIx += 1
            values = np.empty(points.shape)
            for p in range(len(points.flat)):
                values.flat[p] = np.sum(pi0*R**points.flat[p])
            Ret.append(values)
        elif type(argv[argIx]) is str and argv[argIx]=="stDistrPH":
            # transform to ph distribution
            ix = np.arange(N)
            nz = ix[eta.flat>prec]
            Delta = Diag(eta)
            A = np.kron(L+F,I[nz,:][:,nz]) + np.kron(B,Delta[nz,:][:,nz].I*Rh[nz,:][:,nz].T*Delta[nz,:][:,nz])
            alpha = z.T*np.kron(I,Delta[:,nz])
            Ret.append((alpha, A))
        elif type(argv[argIx]) is str and argv[argIx]=="stDistrME":
            # transform it such that the closing vector is a vector of ones
            # this is the way butools accepts ME distributions
            Bm = TransformToOnes(z)
            Bmi = Bm.I
            A = Bm * (np.kron(L.T+F.T,I) + np.kron(B.T,Rh)) * Bmi
            alpha = np.kron(ml.ones((1,N)), eta) * Bmi
            Ret.append((alpha, A))
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
    
    if butools.checkInput and not CheckMAPRepresentation(D0,D1,prec):
        raise Exception('MAPMAP1: The arrival process (D0,D1) is not a valid MAP representation!')
    
    if butools.checkInput and not CheckMAPRepresentation(S0,S1,prec):
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
        elif type(argv[argIx]) is str and argv[argIx]=="qlDistrDPH":
            # transform it to DPH
            alpha = pi0*R*(I-R).I
            A = Diag(alpha).I*R.T*Diag(alpha)
            Ret.append((alpha, A))
        elif type(argv[argIx]) is str and argv[argIx]=="qlDistrMG":
            # transform it to MG
            B = TransformToOnes(np.sum((I-R).I*R,1))
            Bi = B.I
            A = B*R*Bi
            alpha = pi0*Bi
            Ret.append((alpha, A))
        elif type(argv[argIx]) is str and argv[argIx]=="qlMoms":
            numOfMoms = argv[argIx+1]
            argIx += 1
            moms = []
            iR = (I-R).I
            for m in range(1,numOfMoms+1):
                moms.append(math.factorial(m)*np.sum(pi0*iR**(m+1)*R**m))
            Ret.append(MomsFromFactorialMoms(moms))
        elif type(argv[argIx]) is str and argv[argIx]=="qlDistr":
            points = argv[argIx+1]
            argIx += 1
            values = np.empty(points.shape)
            for p in range(len(points.flat)):
                values.flat[p] = np.sum(pi0*R**points.flat[p])
            Ret.append(values)
        elif type(argv[argIx]) is str and argv[argIx]=="stDistrPH":
            # transform it to PH representation
            beta = CTMCSolve(S0+S1, prec)
            theta = DTMCSolve(-D0.I*D1, prec)
            vv = np.kron(theta,beta)
            ix = np.arange(N)
            nz = ix[vv.flat>prec]
            delta = Diag(vv[:,nz])
            alpha = ml.ones((1,N))*B[nz,:].T*delta / np.sum(beta*S1)
            A = delta.I*T[nz,:][:,nz].T*delta
            Ret.append((alpha, A))
        elif type(argv[argIx]) is str and argv[argIx]=="stDistrME":
            Ret.append((eta, T))
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
