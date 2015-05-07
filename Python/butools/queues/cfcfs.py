import numpy as np
import numpy.matlib as ml
import scipy.linalg as la
import butools
import math
from butools.mam import GeneralFluidSolve
from butools.utils import Linsolve, Diag
from butools.reptrans import TransformToOnes
from butools.mc import CTMCSolve, CheckGenerator

def FluidQueue (Q, Rin, Rout, *argv):

    # parse options
    prec = 1e-14
    needST = False
    Q0 = None
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

    if butools.checkInput and not CheckGenerator(Q,False,prec):
        raise Exception('FluidQueue: Generator matrix Q is not Markovian!')

    if butools.checkInput and Q0!=None and not CheckGenerator(Q0,False,prec):
        raise Exception('FluidQueue: Generator matrix Q0 is not Markovian!')

    if butools.checkInput and (np.any(np.diag(Rin)<-prec) or np.any(np.diag(Rout)<-prec)):
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
        elif type(argv[argIx]) is str and argv[argIx]=='qlDistrPH':
            # transform it to PH
            Delta = Diag(Linsolve(K.T,-ini.T)) # Delta = diag (ini*inv(-K));
            A = Delta.I*K.T*Delta
            alpha = np.sum(clo,1).T*Delta
            Ret.append((alpha, A))
        elif type(argv[argIx]) is str and argv[argIx]=='qlDistrME':
            # transform it to ME
            B = TransformToOnes((-K).I*np.sum(clo,1))
            Bi = B.I
            alpha = ini*Bi
            A = B*K*Bi
            Ret.append((alpha, A))
        elif type(argv[argIx]) is str and argv[argIx]=='qlMoms':
            numOfMoms = argv[argIx+1]
            argIx += 1
            moms = []
            iK = -K.I
            for m in range(1,numOfMoms+1):
                moms.append(math.factorial(m)*np.sum(ini*iK**(m+1)*clo))
            Ret.append(moms)
        elif type(argv[argIx]) is str and argv[argIx]=='qlDistr':
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
            Ret.append((alpha, A))
        elif type(argv[argIx]) is str and argv[argIx]=='stDistrME':
            B = TransformToOnes(np.reshape(-K.I*clo*Rin,(N*ini.size,1),'F'))
            Bi = B.I
            alpha = np.kron(ml.ones((1,N)), ini/lambd)*Bi
            A = B*(np.kron(Q.T,ml.eye(K.shape[0])) + np.kron(Rout,K))*Bi        
            Ret.append((alpha, A))
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

    # parse options
    prec = 1e-14
    needST = False
    needQL = False
    Q0 = None
    eaten = []
    for i in range(len(argv)):
        if type(argv[i]) is str and argv[i]=="prec":
            prec = argv[i+1]
            eaten.append(i)
            eaten.append(i+1)
        elif type(argv[i]) is str and len(argv[i])>2 and argv[i][0:2]=="st":
            needST = True
        elif type(argv[i]) is str and len(argv[i])>2 and argv[i][0:2]=="ql":
            needQL = True

    if butools.checkInput and not CheckGenerator(Qin,False,prec):
        raise Exception('FluFluQueue: Generator matrix Qin is not Markovian!')

    if butools.checkInput and not CheckGenerator(Qout,False,prec):
        raise Exception('FluFluQueue: Generator matrix Qout is not Markovian!')

    if butools.checkInput and (np.any(np.diag(Rin)<-prec) or np.any(np.diag(Rout)<-prec)):
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
        massh, inih, Kh, cloh = GeneralFluidSolve (Qh, Rh)

        # sojourn time density in case of 
        # srv0stop = false: inih*expm(Kh*x)*cloh*kron(Rin,Iout)/lambda
        # srv0stop = true: inih*expm(Kh*x)*cloh*kron(Rin,Rout)/lambda/mu    
        lambd = np.sum(CTMCSolve(Qin,prec)*Rin)
        mu = np.sum(CTMCSolve(Qout,prec)*Rout)
    
    Ret = []
    argIx = 0
    while argIx<len(argv):
        if argIx in eaten:
            argIx += 1
            continue
        elif type(argv[argIx]) is str and argv[argIx]=='qlDistrPH':
            # transform it to PH
            Delta = Diag(Linsolve(K.T,-ini.T)) # Delta = diag (ini*inv(-K));
            A = Delta.I*K.T*Delta
            alpha = np.sum(clo,1).T*Delta
            Ret.append((alpha, A))
        elif type(argv[argIx]) is str and argv[argIx]=='qlDistrME':
            # transform it to ME
            B = TransformToOnes((-K).I*np.sum(clo,1))
            Bi = B.I
            alpha = ini*Bi
            A = B*K*Bi
            Ret.append((alpha, A))
        elif type(argv[argIx]) is str and argv[argIx]=='qlMoms':
            numOfMoms = argv[argIx+1]
            argIx += 1
            moms = []
            iK = -K.I
            for m in range(1,numOfMoms+1):
                moms.append(math.factorial(m)*np.sum(ini*iK**(m+1)*clo))
            Ret.append(moms)
        elif type(argv[argIx]) is str and argv[argIx]=='qlDistr':
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
            Ret.append((alpha, A))
        elif type(argv[argIx]) is str and argv[argIx]=='stDistrME':
            # convert result to ME representation
            if not srv0stop:
                B = TransformToOnes(np.sum(cloh*np.kron(Rin,Iout)/lambd,1))
            else:
                B = TransformToOnes(np.sum(cloh*np.kron(Rin,Rout)/lambd/mu,1))
            iB = B.I
            A = B*Kh*iB
            alpha = inih*(-Kh).I*iB
            Ret.append((alpha, A))
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
            values = ml.empty(points.shape)
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
