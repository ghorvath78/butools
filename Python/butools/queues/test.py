import numpy as np
import numpy.matlib as ml
import scipy.linalg as la
import butools
from butools.queues import *
from butools.mc import CheckGenerator, CheckProbMatrix
from butools.dph import *
from butools.ph import *
from butools.map import *

def TestQueuesPackage ():

    import os
    print(os.path.dirname(os.path.realpath(__file__)))

    print("---BuTools: MAM package test file---")
    
    print("Enable the verbose messages with the BuToolsVerbose flag")
    butools.verbose = True
    
    print("Enable input parameter checking with the BuToolsCheckInput flag")
    butools.checkInput = True
    
    # ============================= QBD tests ===============================
    
    print("----------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    # QBD example 1
    B=ml.matrix([[6, 1, 0], [0, 4, 1], [2, 0, 0]])
    F=ml.matrix([[0, 1, 1],[5, 0, 0],[1, 3, 0]])
    L=ml.matrix([[-14, 3, 2], [ 0, -14, 4], [ 3, 1, -10]])
    L0 = L+B
    
    pi0, R = QBDSolve (B, L, F, L0)
    lambd = np.sum(pi0*(ml.eye(R.shape[0])-R).I*F)
    
    print("Test:")
    print("-----")
    
    print("qld,qlm = QBDQueue(B, L, F, L0, 'qlDistr', (0:10), 'qlMoms', 5):")
    qld,qlm = QBDQueue(B, L, F, L0, 'qlDistr', np.arange(0,11), 'qlMoms', 5)
    print("qld=",qld,"\nqlm=",qlm)
    print("std, stm = QBDQueue(B, L, F, L0, 'stDistr', (0:0.1:1), 'stMoms', 5):")
    std, stm = QBDQueue(B, L, F, L0, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5)
    print("std=",std,"\nstm=",stm)
    print("alphap,Ap = QBDQueue(B, L, F, L0, 'qlDistrDPH'):")
    alphap,Ap = QBDQueue(B, L, F, L0, 'qlDistrDPH')
    print("alphap=",alphap,",Ap=",Ap)
    print("alpha,A = QBDQueue(B, L, F, L0, 'qlDistrMG'):")
    alpha,A = QBDQueue(B, L, F, L0, 'qlDistrMG')
    print("alpha=",alpha,",A=",A)
    print("betap,Bp = QBDQueue(B, L, F, L0, 'stDistrPH'):")
    betap, Bp = QBDQueue(B, L, F, L0, 'stDistrPH')
    print("betap=",betap,",Bp=",Bp)
    print("beta,B = QBDQueue(B, L, F, L0, 'stDistrME'):")
    beta, B = QBDQueue(B, L, F, L0, 'stDistrME')
    print("beta=",beta,",B=",B)
     
    assert CheckMGRepresentation(alpha,A), "QBDQueue: invalid MG representation of the queue length!"
    assert CheckMERepresentation(beta,B), "QBDQueue: invalid ME representation of the sojourn time!"
    assert CheckDPHRepresentation(alphap,Ap), "QBDQueue: invalid DPH representation of the queue length!"
    assert CheckPHRepresentation(betap,Bp), "QBDQueue: invalid PH representation of the sojourn time!"
    
    # check Little formula
    mql = MomentsFromMG(alpha,A,1)[0]
    mst = MomentsFromME(beta,B,1)[0]
    assert abs(mql-mst*lambd)<1e-12, "QBDQueue: Little formula does not hold!"
    
    # check the equality of the PH and ME results
    assert la.norm((np.array(MomentsFromDPH(alphap,Ap,5))-np.array(MomentsFromMG(alpha,A,5)))/np.array(MomentsFromMG(alpha,A,5)))<1e-12, "QBDQueue: the MG and DPH representations are not equal!"
    assert la.norm((np.array(MomentsFromPH(betap,Bp,5))-np.array(MomentsFromME(beta,B,5)))/np.array(MomentsFromME(beta,B,5)))<1e-12, "QBDQueue: the ME and PH representations are not equal!"
    
    # check moment and distribution calculation
    assert la.norm(qld-PmfFromMG(alpha,A,np.arange(0,11)))<1e-12, "QBDQueue: qlDistr returns wrong queue length distribution!"
    assert la.norm(std-CdfFromME(beta,B,np.arange(0,1.1,0.1)))<1e-12, "QBDQueue: stDistr returns wrong sojourn time distribution!"
    assert la.norm(np.array(qlm)-np.array(MomentsFromMG(alpha,A,5)))<1e-10, "QBDQueue: qlMoms returns wrong queue length moments!"
    assert la.norm(np.array(stm)-np.array(MomentsFromME(beta,B,5)))<1e-10, "QBDQueue: stMoms returns wrong sojourn time moments!"
    
    print("Input:")
    print("------")
    
    # QBD example 2
    B=ml.matrix([[6, 1, 0], [ 0, 4, 1], [ 2, 0, 0]])
    F=ml.matrix([[0, 0, 0], [ 5, 0, 0], [ 1, 3, 0]])
    L=ml.matrix([[-12, 3, 2], [ 0, -14, 4], [ 3, 1, -10]])
    L0 = L+B

    pi0, R = QBDSolve (B, L, F, L0)
    lambd = np.sum(pi0*(ml.eye(R.shape[0])-R).I*F)
    
    print("Test:")
    print("-----")
    
    print("qld,qlm = QBDQueue(B, L, F, L0, 'qlDistr', (0:10), 'qlMoms', 5):")
    qld,qlm = QBDQueue(B, L, F, L0, 'qlDistr', np.arange(0,11), 'qlMoms', 5)
    print("qld=",qld,"\nqlm=",qlm)
    print("std, stm = QBDQueue(B, L, F, L0, 'stDistr', (0:0.1:1), 'stMoms', 5):")
    std, stm = QBDQueue(B, L, F, L0, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5)
    print("std=",std,"\nstm=",stm)
    print("alphap,Ap = QBDQueue(B, L, F, L0, 'qlDistrDPH'):")
    alphap,Ap = QBDQueue(B, L, F, L0, 'qlDistrDPH')
    print("alphap=",alphap,",Ap=",Ap)
    print("alpha,A = QBDQueue(B, L, F, L0, 'qlDistrMG'):")
    alpha,A = QBDQueue(B, L, F, L0, 'qlDistrMG')
    print("alpha=",alpha,",A=",A)
    print("betap,Bp = QBDQueue(B, L, F, L0, 'stDistrPH'):")
    betap, Bp = QBDQueue(B, L, F, L0, 'stDistrPH')
    print("betap=",betap,",Bp=",Bp)
    print("beta,B = QBDQueue(B, L, F, L0, 'stDistrME'):")
    beta, B = QBDQueue(B, L, F, L0, 'stDistrME')
    print("beta=",beta,",B=",B)
     
    assert CheckMGRepresentation(alpha,A), "QBDQueue: invalid MG representation of the queue length!"
    assert CheckMERepresentation(beta,B), "QBDQueue: invalid ME representation of the sojourn time!"
    assert CheckDPHRepresentation(alphap,Ap), "QBDQueue: invalid DPH representation of the queue length!"
    assert CheckPHRepresentation(betap,Bp), "QBDQueue: invalid PH representation of the sojourn time!"
    
    # check Little formula
    mql = MomentsFromMG(alpha,A,1)[0]
    mst = MomentsFromME(beta,B,1)[0]
    assert abs(mql-mst*lambd)<1e-12, "QBDQueue: Little formula does not hold!"
    
    # check the equality of the PH and ME results
    assert la.norm((np.array(MomentsFromDPH(alphap,Ap,5))-np.array(MomentsFromMG(alpha,A,5)))/np.array(MomentsFromMG(alpha,A,5)))<1e-12, "QBDQueue: the MG and DPH representations are not equal!"
    assert la.norm((np.array(MomentsFromPH(betap,Bp,5))-np.array(MomentsFromME(beta,B,5)))/np.array(MomentsFromME(beta,B,5)))<1e-12, "QBDQueue: the ME and PH representations are not equal!"
    
    # check moment and distribution calculation
    assert la.norm(qld-PmfFromMG(alpha,A,np.arange(0,11)))<1e-12, "QBDQueue: qlDistr returns wrong queue length distribution!"
    assert la.norm(std-CdfFromME(beta,B,np.arange(0,1.1,0.1)))<1e-12, "QBDQueue: stDistr returns wrong sojourn time distribution!"
    assert la.norm(np.array(qlm)-np.array(MomentsFromMG(alpha,A,5)))<1e-10, "QBDQueue: qlMoms returns wrong queue length moments!"
    assert la.norm(np.array(stm)-np.array(MomentsFromME(beta,B,5)))<1e-10, "QBDQueue: stMoms returns wrong sojourn time moments!"
    
    print("Input:")
    print("------")
    
    # QBD example 3
    B=ml.matrix([[6, 1, 0], [ 0, 5, 0], [ 0, 0, 0]])
    F=ml.matrix([[0, 3, 1], [ 0, 5, 0], [ 0, 0, 0]])
    L=ml.matrix([[-16, 3, 2], [ 0, -14, 4], [ 3, 1, -4]])
    L0=ml.matrix([[-14, 10, 0], [ 5, -10, 0], [ 3, 1, -4]])
    
    pi0, R = QBDSolve (B, L, F, L0)
    lambd = np.sum(pi0*(ml.eye(R.shape[0])-R).I*F)
    
    print("Test:")
    print("-----")
    
    print("qld,qlm = QBDQueue(B, L, F, L0, 'qlDistr', (0:10), 'qlMoms', 5):")
    qld,qlm = QBDQueue(B, L, F, L0, 'qlDistr', np.arange(0,11), 'qlMoms', 5)
    print("qld=",qld,"\nqlm=",qlm)
    print("std, stm = QBDQueue(B, L, F, L0, 'stDistr', (0:0.1:1), 'stMoms', 5):")
    std, stm = QBDQueue(B, L, F, L0, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5)
    print("std=",std,"\nstm=",stm)
    print("alphap,Ap = QBDQueue(B, L, F, L0, 'qlDistrDPH'):")
    alphap,Ap = QBDQueue(B, L, F, L0, 'qlDistrDPH')
    print("alphap=",alphap,",Ap=",Ap)
    print("alpha,A = QBDQueue(B, L, F, L0, 'qlDistrMG'):")
    alpha,A = QBDQueue(B, L, F, L0, 'qlDistrMG')
    print("alpha=",alpha,",A=",A)
    print("betap,Bp = QBDQueue(B, L, F, L0, 'stDistrPH'):")
    betap, Bp = QBDQueue(B, L, F, L0, 'stDistrPH')
    print("betap=",betap,",Bp=",Bp)
    print("beta,B = QBDQueue(B, L, F, L0, 'stDistrME'):")
    beta, B = QBDQueue(B, L, F, L0, 'stDistrME')
    print("beta=",beta,",B=",B)
     
    assert CheckMGRepresentation(alpha,A), "QBDQueue: invalid MG representation of the queue length!"
    assert CheckMERepresentation(beta,B), "QBDQueue: invalid ME representation of the sojourn time!"
    assert CheckDPHRepresentation(alphap,Ap), "QBDQueue: invalid DPH representation of the queue length!"
    assert CheckPHRepresentation(betap,Bp), "QBDQueue: invalid PH representation of the sojourn time!"
    
    # check Little formula
    mql = MomentsFromMG(alpha,A,1)[0]
    mst = MomentsFromME(beta,B,1)[0]
    assert abs(mql-mst*lambd)<1e-12, "QBDQueue: Little formula does not hold!"
    
    # check the equality of the PH and ME results
    assert la.norm((np.array(MomentsFromDPH(alphap,Ap,5))-np.array(MomentsFromMG(alpha,A,5)))/np.array(MomentsFromMG(alpha,A,5)))<1e-12, "QBDQueue: the MG and DPH representations are not equal!"
    assert la.norm((np.array(MomentsFromPH(betap,Bp,5))-np.array(MomentsFromME(beta,B,5)))/np.array(MomentsFromME(beta,B,5)))<1e-12, "QBDQueue: the ME and PH representations are not equal!"
    
    # check moment and distribution calculation
    assert la.norm(qld-PmfFromMG(alpha,A,np.arange(0,11)))<1e-12, "QBDQueue: qlDistr returns wrong queue length distribution!"
    assert la.norm(std-CdfFromME(beta,B,np.arange(0,1.1,0.1)))<1e-12, "QBDQueue: stDistr returns wrong sojourn time distribution!"
    assert la.norm(np.array(qlm)-np.array(MomentsFromMG(alpha,A,5)))<1e-7, "QBDQueue: qlMoms returns wrong queue length moments!"
    assert la.norm(np.array(stm)-np.array(MomentsFromME(beta,B,5)))<1e-7, "QBDQueue: stMoms returns wrong sojourn time moments!"
    
    print("Input:")
    print("------")
    
    # QBD example 4
    B = ml.matrix([[0,0], [ 3,4]])
    L = ml.matrix([[-6,5], [ 3,-12]])
    F = ml.matrix([[1,0], [ 2,0]])
    L0 = ml.matrix([[-6,5], [ 6,-8]])

    pi0, R = QBDSolve (B, L, F, L0)
    lambd = np.sum(pi0*(ml.eye(R.shape[0])-R).I*F)
    
    print("Test:")
    print("-----")
    
    print("qld,qlm = QBDQueue(B, L, F, L0, 'qlDistr', (0:10), 'qlMoms', 5):")
    qld,qlm = QBDQueue(B, L, F, L0, 'qlDistr', np.arange(0,11), 'qlMoms', 5)
    print("qld=",qld,"\nqlm=",qlm)
    print("std, stm = QBDQueue(B, L, F, L0, 'stDistr', (0:0.1:1), 'stMoms', 5):")
    std, stm = QBDQueue(B, L, F, L0, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5)
    print("std=",std,"\nstm=",stm)
    print("alphap,Ap = QBDQueue(B, L, F, L0, 'qlDistrDPH'):")
    alphap,Ap = QBDQueue(B, L, F, L0, 'qlDistrDPH')
    print("alphap=",alphap,",Ap=",Ap)
    print("alpha,A = QBDQueue(B, L, F, L0, 'qlDistrMG'):")
    alpha,A = QBDQueue(B, L, F, L0, 'qlDistrMG')
    print("alpha=",alpha,",A=",A)
    print("betap,Bp = QBDQueue(B, L, F, L0, 'stDistrPH'):")
    betap, Bp = QBDQueue(B, L, F, L0, 'stDistrPH')
    print("betap=",betap,",Bp=",Bp)
    print("beta,B = QBDQueue(B, L, F, L0, 'stDistrME'):")
    beta, B = QBDQueue(B, L, F, L0, 'stDistrME')
    print("beta=",beta,",B=",B)
     
    assert CheckMGRepresentation(alpha,A), "QBDQueue: invalid MG representation of the queue length!"
    assert CheckMERepresentation(beta,B), "QBDQueue: invalid ME representation of the sojourn time!"
    assert CheckDPHRepresentation(alphap,Ap), "QBDQueue: invalid DPH representation of the queue length!"
    assert CheckPHRepresentation(betap,Bp), "QBDQueue: invalid PH representation of the sojourn time!"
    
    # check Little formula
    mql = MomentsFromMG(alpha,A,1)[0]
    mst = MomentsFromME(beta,B,1)[0]
    assert abs(mql-mst*lambd)<1e-12, "QBDQueue: Little formula does not hold!"
    
    # check the equality of the PH and ME results
    assert la.norm((np.array(MomentsFromDPH(alphap,Ap,5))-np.array(MomentsFromMG(alpha,A,5)))/np.array(MomentsFromMG(alpha,A,5)))<1e-12, "QBDQueue: the MG and DPH representations are not equal!"
    assert la.norm((np.array(MomentsFromPH(betap,Bp,5))-np.array(MomentsFromME(beta,B,5)))/np.array(MomentsFromME(beta,B,5)))<1e-12, "QBDQueue: the ME and PH representations are not equal!"
    
    # check moment and distribution calculation
    assert la.norm(qld-PmfFromMG(alpha,A,np.arange(0,11)))<1e-12, "QBDQueue: qlDistr returns wrong queue length distribution!"
    assert la.norm(std-CdfFromME(beta,B,np.arange(0,1.1,0.1)))<1e-12, "QBDQueue: stDistr returns wrong sojourn time distribution!"
    assert la.norm(np.array(qlm)-np.array(MomentsFromMG(alpha,A,5)))<1e-7, "QBDQueue: qlMoms returns wrong queue length moments!"
    assert la.norm(np.array(stm)-np.array(MomentsFromME(beta,B,5)))<1e-7, "QBDQueue: stMoms returns wrong sojourn time moments!"

    # ============================= MAP/MAP/1 tests ===============================
    
    print("----------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-8, 2], [ 1, -3]])
    D1=ml.matrix([[1, 5], [ 0, 2]])
    
    S0=ml.matrix([[-10, 4], [ 0, -7]])
    S1=ml.matrix([[5, 1], [ 4, 3]])
    
    print("Test:")
    print("-----")

    print("qld,qlm = MAPMAP1(D0,D1,S0,S1, 'qlDistr', (0:10), 'qlMoms', 5):")
    qld,qlm = MAPMAP1(D0,D1,S0,S1, 'qlDistr', np.arange(0,11), 'qlMoms', 5)
    print("qld=",qld,"\nqlm=",qlm)
    print("std, stm = MAPMAP1(D0,D1,S0,S1, 'stDistr', (0:0.1:1), 'stMoms', 5):")
    std, stm = MAPMAP1(D0,D1,S0,S1, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5)
    print("std=",std,"\nstm=",stm)
    print("alphap,Ap = MAPMAP1(D0,D1,S0,S1, 'qlDistrDPH'):")
    alphap,Ap = MAPMAP1(D0,D1,S0,S1, 'qlDistrDPH')
    print("alphap=",alphap,",Ap=",Ap)
    print("alpha,A = MAPMAP1(D0,D1,S0,S1, 'qlDistrMG'):")
    alpha,A = MAPMAP1(D0,D1,S0,S1, 'qlDistrMG')
    print("alpha=",alpha,",A=",A)
    print("betap,Bp = MAPMAP1(D0,D1,S0,S1, 'stDistrPH'):")
    betap, Bp = MAPMAP1(D0,D1,S0,S1, 'stDistrPH')
    print("betap=",betap,",Bp=",Bp)
    print("beta,B = MAPMAP1(D0,D1,S0,S1, 'stDistrME'):")
    beta, B = MAPMAP1(D0,D1,S0,S1, 'stDistrME')
    print("beta=",beta,",B=",B)
     
    assert CheckMGRepresentation(alpha,A), "MAPMAP1: invalid MG representation of the queue length!"
    assert CheckMERepresentation(beta,B), "MAPMAP1: invalid ME representation of the sojourn time!"
    assert CheckDPHRepresentation(alphap,Ap), "MAPMAP1: invalid DPH representation of the queue length!"
    assert CheckPHRepresentation(betap,Bp), "MAPMAP1: invalid PH representation of the sojourn time!"
    
    # cross-check
    IA = ml.eye(D0.shape[0])
    IS = ml.eye(S0.shape[0])
    gamma, G = QBDQueue (np.kron(IA,S1), np.kron(D0,IS)+np.kron(IA,S0), np.kron(D1,IS), np.kron(D0,IS), "stDistrME")
    msmall = np.array(MomentsFromME(beta,B,5))
    mlarge = np.array(MomentsFromME(gamma,G,5))
    assert la.norm((msmall-mlarge)/msmall)<1e-12, "MAPMAP1: Large and small model does not give the same results!"

    # check Little formula
    mql = MomentsFromMG(alpha,A,1)[0]
    mst = MomentsFromME(beta,B,1)[0]
    lambd = 1/MarginalMomentsFromMAP(D0,D1,1)[0]
    assert abs(mql-mst*lambd)<1e-12, "MAPMAP1: Little formula does not hold!"
    
    # check the equality of the PH and ME results
    assert la.norm((np.array(MomentsFromDPH(alphap,Ap,5))-np.array(MomentsFromMG(alpha,A,5)))/np.array(MomentsFromMG(alpha,A,5)))<1e-12, "MAPMAP1: the MG and DPH representations are not equal!"
    assert la.norm((np.array(MomentsFromPH(betap,Bp,5))-np.array(MomentsFromME(beta,B,5)))/np.array(MomentsFromME(beta,B,5)))<1e-12, "MAPMAP1: the ME and PH representations are not equal!"
    
    # check moment and distribution calculation
    assert la.norm(qld-PmfFromMG(alpha,A,np.arange(0,11)))<1e-12, "MAPMAP1: qlDistr returns wrong queue length distribution!"
    assert la.norm(std-CdfFromME(beta,B,np.arange(0,1.1,0.1)))<1e-12, "MAPMAP1: stDistr returns wrong sojourn time distribution!"
    assert la.norm(np.array(qlm)-np.array(MomentsFromMG(alpha,A,5)))<1e-7, "MAPMAP1: qlMoms returns wrong queue length moments!"
    assert la.norm(np.array(stm)-np.array(MomentsFromME(beta,B,5)))<1e-7, "MAPMAP1: stMoms returns wrong sojourn time moments!"
    
    
    print("Input:")
    print("------")
    
    D0 = ml.matrix([[-8, 1, 2], [ 0, -6, 4], [ 3, 0, -3]])
    D1 = ml.matrix([[4, 1, 0], [ 0, 2, 0], [ 0, 0, 0]])

    print("Test:")
    print("-----")

    print("qld,qlm = MAPMAP1(D0,D1,S0,S1, 'qlDistr', (0:10), 'qlMoms', 5):")
    qld,qlm = MAPMAP1(D0,D1,S0,S1, 'qlDistr', np.arange(0,11), 'qlMoms', 5)
    print("qld=",qld,"\nqlm=",qlm)
    print("std, stm = MAPMAP1(D0,D1,S0,S1, 'stDistr', (0:0.1:1), 'stMoms', 5):")
    std, stm = MAPMAP1(D0,D1,S0,S1, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5)
    print("std=",std,"\nstm=",stm)
    print("alphap,Ap = MAPMAP1(D0,D1,S0,S1, 'qlDistrDPH'):")
    alphap,Ap = MAPMAP1(D0,D1,S0,S1, 'qlDistrDPH')
    print("alphap=",alphap,",Ap=",Ap)
    print("alpha,A = MAPMAP1(D0,D1,S0,S1, 'qlDistrMG'):")
    alpha,A = MAPMAP1(D0,D1,S0,S1, 'qlDistrMG')
    print("alpha=",alpha,",A=",A)
    print("betap,Bp = MAPMAP1(D0,D1,S0,S1, 'stDistrPH'):")
    betap, Bp = MAPMAP1(D0,D1,S0,S1, 'stDistrPH')
    print("betap=",betap,",Bp=",Bp)
    print("beta,B = MAPMAP1(D0,D1,S0,S1, 'stDistrME'):")
    beta, B = MAPMAP1(D0,D1,S0,S1, 'stDistrME')
    print("beta=",beta,",B=",B)
     
    assert CheckMGRepresentation(alpha,A), "MAPMAP1: invalid MG representation of the queue length!"
    assert CheckMERepresentation(beta,B), "MAPMAP1: invalid ME representation of the sojourn time!"
    assert CheckDPHRepresentation(alphap,Ap), "MAPMAP1: invalid DPH representation of the queue length!"
    assert CheckPHRepresentation(betap,Bp), "MAPMAP1: invalid PH representation of the sojourn time!"
    
    # cross-check
    IA = ml.eye(D0.shape[0])
    IS = ml.eye(S0.shape[0])
    gamma, G = QBDQueue (np.kron(IA,S1), np.kron(D0,IS)+np.kron(IA,S0), np.kron(D1,IS), np.kron(D0,IS), "stDistrME")
    msmall = np.array(MomentsFromME(beta,B,5))
    mlarge = np.array(MomentsFromME(gamma,G,5))
    assert la.norm((msmall-mlarge)/msmall)<1e-12, "MAPMAP1: Large and small model does not give the same results!"

    # check Little formula
    mql = MomentsFromMG(alpha,A,1)[0]
    mst = MomentsFromME(beta,B,1)[0]
    lambd = 1/MarginalMomentsFromMAP(D0,D1,1)[0]
    assert abs(mql-mst*lambd)<1e-12, "MAPMAP1: Little formula does not hold!"
    
    # check the equality of the PH and ME results
    assert la.norm((np.array(MomentsFromDPH(alphap,Ap,5))-np.array(MomentsFromMG(alpha,A,5)))/np.array(MomentsFromMG(alpha,A,5)))<1e-12, "MAPMAP1: the MG and DPH representations are not equal!"
    assert la.norm((np.array(MomentsFromPH(betap,Bp,5))-np.array(MomentsFromME(beta,B,5)))/np.array(MomentsFromME(beta,B,5)))<1e-12, "MAPMAP1: the ME and PH representations are not equal!"
    
    # check moment and distribution calculation
    assert la.norm(qld-PmfFromMG(alpha,A,np.arange(0,11)))<1e-12, "MAPMAP1: qlDistr returns wrong queue length distribution!"
    assert la.norm(std-CdfFromME(beta,B,np.arange(0,1.1,0.1)))<1e-12, "MAPMAP1: stDistr returns wrong sojourn time distribution!"
    assert la.norm(np.array(qlm)-np.array(MomentsFromMG(alpha,A,5)))<1e-7, "MAPMAP1: qlMoms returns wrong queue length moments!"
    assert la.norm(np.array(stm)-np.array(MomentsFromME(beta,B,5)))<1e-7, "MAPMAP1: stMoms returns wrong sojourn time moments!"

    
    print("Input:")
    print("------")
    
    S0 = ml.matrix([[-10, 4, 0], [ 5, -7, 2], [ 1, 2, -8]])
    S1 = ml.matrix([[0, 0, 6], [ 0, 0, 0], [ 0, 3, 2]])
    
    D0 = ml.matrix([[-8, 1, 2], [ 0, -6, 4], [ 3, 0, -3]])
    D1 = ml.matrix([[4, 1, 0], [ 0, 0, 2], [ 0, 0, 0]])

    print("Test:")
    print("-----")

    print("qld,qlm = MAPMAP1(D0,D1,S0,S1, 'qlDistr', (0:10), 'qlMoms', 5):")
    qld,qlm = MAPMAP1(D0,D1,S0,S1, 'qlDistr', np.arange(0,11), 'qlMoms', 5)
    print("qld=",qld,"\nqlm=",qlm)
    print("std, stm = MAPMAP1(D0,D1,S0,S1, 'stDistr', (0:0.1:1), 'stMoms', 5):")
    std, stm = MAPMAP1(D0,D1,S0,S1, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5)
    print("std=",std,"\nstm=",stm)
    print("alphap,Ap = MAPMAP1(D0,D1,S0,S1, 'qlDistrDPH'):")
    alphap,Ap = MAPMAP1(D0,D1,S0,S1, 'qlDistrDPH')
    print("alphap=",alphap,",Ap=",Ap)
    print("alpha,A = MAPMAP1(D0,D1,S0,S1, 'qlDistrMG'):")
    alpha,A = MAPMAP1(D0,D1,S0,S1, 'qlDistrMG')
    print("alpha=",alpha,",A=",A)
    print("betap,Bp = MAPMAP1(D0,D1,S0,S1, 'stDistrPH'):")
    betap, Bp = MAPMAP1(D0,D1,S0,S1, 'stDistrPH')
    print("betap=",betap,",Bp=",Bp)
    print("beta,B = MAPMAP1(D0,D1,S0,S1, 'stDistrME'):")
    beta, B = MAPMAP1(D0,D1,S0,S1, 'stDistrME')
    print("beta=",beta,",B=",B)
     
    assert CheckMGRepresentation(alpha,A), "MAPMAP1: invalid MG representation of the queue length!"
    assert CheckMERepresentation(beta,B), "MAPMAP1: invalid ME representation of the sojourn time!"
    assert CheckDPHRepresentation(alphap,Ap), "MAPMAP1: invalid DPH representation of the queue length!"
    assert CheckPHRepresentation(betap,Bp), "MAPMAP1: invalid PH representation of the sojourn time!"
    
    # cross-check
    IA = ml.eye(D0.shape[0])
    IS = ml.eye(S0.shape[0])
    gamma, G = QBDQueue (np.kron(IA,S1), np.kron(D0,IS)+np.kron(IA,S0), np.kron(D1,IS), np.kron(D0,IS), "stDistrME")
    msmall = np.array(MomentsFromME(beta,B,5))
    mlarge = np.array(MomentsFromME(gamma,G,5))
    assert la.norm((msmall-mlarge)/msmall)<1e-12, "MAPMAP1: Large and small model does not give the same results!"

    # check Little formula
    mql = MomentsFromMG(alpha,A,1)[0]
    mst = MomentsFromME(beta,B,1)[0]
    lambd = 1/MarginalMomentsFromMAP(D0,D1,1)[0]
    assert abs(mql-mst*lambd)<1e-12, "MAPMAP1: Little formula does not hold!"
    
    # check the equality of the PH and ME results
    assert la.norm((np.array(MomentsFromDPH(alphap,Ap,5))-np.array(MomentsFromMG(alpha,A,5)))/np.array(MomentsFromMG(alpha,A,5)))<1e-12, "MAPMAP1: the MG and DPH representations are not equal!"
    assert la.norm((np.array(MomentsFromPH(betap,Bp,5))-np.array(MomentsFromME(beta,B,5)))/np.array(MomentsFromME(beta,B,5)))<1e-12, "MAPMAP1: the ME and PH representations are not equal!"
    
    # check moment and distribution calculation
    assert la.norm(qld-PmfFromMG(alpha,A,np.arange(0,11)))<1e-12, "MAPMAP1: qlDistr returns wrong queue length distribution!"
    assert la.norm(std-CdfFromME(beta,B,np.arange(0,1.1,0.1)))<1e-12, "MAPMAP1: stDistr returns wrong sojourn time distribution!"
    assert la.norm(np.array(qlm)-np.array(MomentsFromMG(alpha,A,5)))<1e-7, "MAPMAP1: qlMoms returns wrong queue length moments!"
    assert la.norm(np.array(stm)-np.array(MomentsFromME(beta,B,5)))<1e-7, "MAPMAP1: stMoms returns wrong sojourn time moments!"

    
    # ============================= MAP/PH/1 tests ===============================
    
    print("----------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0 = ml.matrix([[-8, 1, 2], [ 0, -6, 4], [ 3, 0, -3]])
    D1 = ml.matrix([[4, 1, 0], [ 0, 0, 2], [ 0, 0, 0]])
    sigma = ml.matrix([[0.2, 0.7, 0.1]])
    S = ml.matrix([[-10, 4, 0], [ 5, -7, 2], [ 1, 2, -8]])
    
    S0 = S
    S1 = np.sum(-S,1)*sigma

    print("Test:")
    print("-----")

    print("qld,qlm = MAPMAP1(D0,D1,S0,S1, 'qlDistr', (0:10), 'qlMoms', 5):")
    qld,qlm = MAPMAP1(D0,D1,S0,S1, 'qlDistr', np.arange(0,11), 'qlMoms', 5)
    print("qld=",qld,"\nqlm=",qlm)
    print("std, stm = MAPMAP1(D0,D1,S0,S1, 'stDistr', (0:0.1:1), 'stMoms', 5):")
    std, stm = MAPMAP1(D0,D1,S0,S1, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5)
    print("std=",std,"\nstm=",stm)
    print("alphap,Ap = MAPMAP1(D0,D1,S0,S1, 'qlDistrDPH'):")
    alphap,Ap = MAPMAP1(D0,D1,S0,S1, 'qlDistrDPH')
    print("alphap=",alphap,",Ap=",Ap)
    print("alpha,A = MAPMAP1(D0,D1,S0,S1, 'qlDistrMG'):")
    alpha,A = MAPMAP1(D0,D1,S0,S1, 'qlDistrMG')
    print("alpha=",alpha,",A=",A)
    print("betap,Bp = MAPMAP1(D0,D1,S0,S1, 'stDistrPH'):")
    betap, Bp = MAPMAP1(D0,D1,S0,S1, 'stDistrPH')
    print("betap=",betap,",Bp=",Bp)
    print("beta,B = MAPMAP1(D0,D1,S0,S1, 'stDistrME'):")
    beta, B = MAPMAP1(D0,D1,S0,S1, 'stDistrME')
    print("beta=",beta,",B=",B)
     
    assert CheckMGRepresentation(alpha,A), "MAPMAP1: invalid MG representation of the queue length!"
    assert CheckMERepresentation(beta,B), "MAPMAP1: invalid ME representation of the sojourn time!"
    assert CheckDPHRepresentation(alphap,Ap), "MAPMAP1: invalid DPH representation of the queue length!"
    assert CheckPHRepresentation(betap,Bp), "MAPMAP1: invalid PH representation of the sojourn time!"
    
    # cross-check
    IA = ml.eye(D0.shape[0])
    IS = ml.eye(S0.shape[0])
    gamma, G = QBDQueue (np.kron(IA,S1), np.kron(D0,IS)+np.kron(IA,S0), np.kron(D1,IS), np.kron(D0,IS), "stDistrME")
    msmall = np.array(MomentsFromME(beta,B,5))
    mlarge = np.array(MomentsFromME(gamma,G,5))
    assert la.norm((msmall-mlarge)/msmall)<1e-12, "MAPMAP1: Large and small model does not give the same results!"

    # check Little formula
    mql = MomentsFromMG(alpha,A,1)[0]
    mst = MomentsFromME(beta,B,1)[0]
    lambd = 1/MarginalMomentsFromMAP(D0,D1,1)[0]
    assert abs(mql-mst*lambd)<1e-12, "MAPMAP1: Little formula does not hold!"
    
    # check the equality of the PH and ME results
    assert la.norm((np.array(MomentsFromDPH(alphap,Ap,5))-np.array(MomentsFromMG(alpha,A,5)))/np.array(MomentsFromMG(alpha,A,5)))<1e-12, "MAPMAP1: the MG and DPH representations are not equal!"
    assert la.norm((np.array(MomentsFromPH(betap,Bp,5))-np.array(MomentsFromME(beta,B,5)))/np.array(MomentsFromME(beta,B,5)))<1e-12, "MAPMAP1: the ME and PH representations are not equal!"
    
    # check moment and distribution calculation
    assert la.norm(qld-PmfFromMG(alpha,A,np.arange(0,11)))<1e-12, "MAPMAP1: qlDistr returns wrong queue length distribution!"
    assert la.norm(std-CdfFromME(beta,B,np.arange(0,1.1,0.1)))<1e-12, "MAPMAP1: stDistr returns wrong sojourn time distribution!"
    assert la.norm(np.array(qlm)-np.array(MomentsFromMG(alpha,A,5)))<1e-7, "MAPMAP1: qlMoms returns wrong queue length moments!"
    assert la.norm(np.array(stm)-np.array(MomentsFromME(beta,B,5)))<1e-7, "MAPMAP1: stMoms returns wrong sojourn time moments!"


    # ============================= PH/PH/1 tests ===============================
    
    print("----------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    delta = ml.matrix([[0.5, 0.1, 0.4]])
    D = ml.matrix([[-8, 1, 2], [ 0, -6, 4], [ 3, 0, -3]])
    sigma = ml.matrix([[0.2, 0.7, 0.1]])
    S = ml.matrix([[-10, 4, 0], [ 5, -7, 2], [ 1, 2, -8]])
    
    D0 = D
    D1 = np.sum(-D,1)*delta
    S0 = S
    S1 = np.sum(-S,1)*sigma

    print("Test:")
    print("-----")

    print("qld,qlm = MAPMAP1(D0,D1,S0,S1, 'qlDistr', (0:10), 'qlMoms', 5):")
    qld,qlm = MAPMAP1(D0,D1,S0,S1, 'qlDistr', np.arange(0,11), 'qlMoms', 5)
    print("qld=",qld,"\nqlm=",qlm)
    print("std, stm = MAPMAP1(D0,D1,S0,S1, 'stDistr', (0:0.1:1), 'stMoms', 5):")
    std, stm = MAPMAP1(D0,D1,S0,S1, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5)
    print("std=",std,"\nstm=",stm)
    print("alphap,Ap = MAPMAP1(D0,D1,S0,S1, 'qlDistrDPH'):")
    alphap,Ap = MAPMAP1(D0,D1,S0,S1, 'qlDistrDPH')
    print("alphap=",alphap,",Ap=",Ap)
    print("alpha,A = MAPMAP1(D0,D1,S0,S1, 'qlDistrMG'):")
    alpha,A = MAPMAP1(D0,D1,S0,S1, 'qlDistrMG')
    print("alpha=",alpha,",A=",A)
    print("betap,Bp = MAPMAP1(D0,D1,S0,S1, 'stDistrPH'):")
    betap, Bp = MAPMAP1(D0,D1,S0,S1, 'stDistrPH')
    print("betap=",betap,",Bp=",Bp)
    print("beta,B = MAPMAP1(D0,D1,S0,S1, 'stDistrME'):")
    beta, B = MAPMAP1(D0,D1,S0,S1, 'stDistrME')
    print("beta=",beta,",B=",B)
     
    assert CheckMGRepresentation(alpha,A), "MAPMAP1: invalid MG representation of the queue length!"
    assert CheckMERepresentation(beta,B), "MAPMAP1: invalid ME representation of the sojourn time!"
    assert CheckDPHRepresentation(alphap,Ap), "MAPMAP1: invalid DPH representation of the queue length!"
    assert CheckPHRepresentation(betap,Bp), "MAPMAP1: invalid PH representation of the sojourn time!"
    
    # cross-check
    IA = ml.eye(D0.shape[0])
    IS = ml.eye(S0.shape[0])
    gamma, G = QBDQueue (np.kron(IA,S1), np.kron(D0,IS)+np.kron(IA,S0), np.kron(D1,IS), np.kron(D0,IS), "stDistrME")
    msmall = np.array(MomentsFromME(beta,B,5))
    mlarge = np.array(MomentsFromME(gamma,G,5))
    assert la.norm((msmall-mlarge)/msmall)<1e-12, "MAPMAP1: Large and small model does not give the same results!"

    # check Little formula
    mql = MomentsFromMG(alpha,A,1)[0]
    mst = MomentsFromME(beta,B,1)[0]
    lambd = 1/MarginalMomentsFromMAP(D0,D1,1)[0]
    assert abs(mql-mst*lambd)<1e-12, "MAPMAP1: Little formula does not hold!"
    
    # check the equality of the PH and ME results
    assert la.norm((np.array(MomentsFromDPH(alphap,Ap,5))-np.array(MomentsFromMG(alpha,A,5)))/np.array(MomentsFromMG(alpha,A,5)))<1e-12, "MAPMAP1: the MG and DPH representations are not equal!"
    assert la.norm((np.array(MomentsFromPH(betap,Bp,5))-np.array(MomentsFromME(beta,B,5)))/np.array(MomentsFromME(beta,B,5)))<1e-12, "MAPMAP1: the ME and PH representations are not equal!"
    
    # check moment and distribution calculation
    assert la.norm(qld-PmfFromMG(alpha,A,np.arange(0,11)))<1e-12, "MAPMAP1: qlDistr returns wrong queue length distribution!"
    assert la.norm(std-CdfFromME(beta,B,np.arange(0,1.1,0.1)))<1e-12, "MAPMAP1: stDistr returns wrong sojourn time distribution!"
    assert la.norm(np.array(qlm)-np.array(MomentsFromMG(alpha,A,5)))<1e-7, "MAPMAP1: qlMoms returns wrong queue length moments!"
    assert la.norm(np.array(stm)-np.array(MomentsFromME(beta,B,5)))<1e-7, "MAPMAP1: stMoms returns wrong sojourn time moments!"

    
    # ============================= MMAPml.matrix([[K]])/PHml.matrix([[K]])/1 PR tests ============================
    
    print("----------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    D0=ml.matrix([[-5.49, 0, 1.15, 0], [ 0, -2.29, 0, 0], [ 0, 0.08, -1.32, 0], [ 0.72, 1.17, 0.7, -7.07]])
    D1=ml.matrix([[0.25, 0.38, 0.64, 0], [ 0, 0, 0, 1.09], [ 0, 1.24, 0, 0], [ 0.37, 0, 0, 0]])
    D2=ml.matrix([[0.3, 1.0, 0, 0.48], [ 0, 0.2, 0, 0], [ 0, 0, 0, 0], [ 0.61, 0, 0, 0.2]])
    D3=ml.matrix([[0, 0.98, 0, 0.31], [ 0, 0, 1.0, 0], [ 0, 0, 0, 0], [ 1.1, 0.84, 0.33, 1.03]])
    
    sigma3 = ml.matrix([[1/6, 5/6]])
    S3 = ml.matrix([[-0.5, 0.5], [ 0, -3]])
    sigma2 = ml.matrix([[1/1.7, 1-1/1.7]])
    S2 = ml.matrix([[-2/0.85, 2/0.85], [ 0, -4]])
    sigma1 = ml.matrix([[1/4, 1-1/4]])
    S1 = ml.matrix([[-2.5, 2.5], [ 0, -10]])
    
    print("Test:")
    print("-----")
    
    qlm, qld = MMAPPH1PRPR((D0,D1,D2,D3), (sigma1,sigma2,sigma3),(S1,S2,S3), 'qlMoms', 3, 'qlDistr', 500)
    momFromDistr = ml.vstack((ml.matrix(np.arange(500))*qld,ml.matrix(np.arange(500)**2)*qld,ml.matrix(np.arange(500)**3)*qld))
    assert la.norm((momFromDistr-qlm)/qlm)<0.001, "MMAPPH1PRPR: queue length moments and queue length distribution are not consistent!"
    
    stm, std = MMAPPH1PRPR((D0,D1,D2,D3), (sigma1,sigma2,sigma3),(S1,S2,S3), 'stMoms', 3, 'stDistr', [1,5,10])   
    assert np.min(std)>=0 and np.max(std)<=1 and np.all(np.diff(std)>=0), "MMAPPH1PRPR: invalid sojourn time distribution!"

    # check Little formula for every traffic class
    lambda1 = 1/MarginalMomentsFromMAP(D0+D2+D3,D1,1)[0]
    lambda2 = 1/MarginalMomentsFromMAP(D0+D1+D3,D2,1)[0]
    lambda3 = 1/MarginalMomentsFromMAP(D0+D1+D2,D3,1)[0]
    assert abs(qlm[0,0]-stm[0,0]*lambda1)<1e-12, "MMAPPH1PRPR: Little formula does not hold for class 1!"
    assert abs(qlm[0,1]-stm[0,1]*lambda2)<1e-12, "MMAPPH1PRPR: Little formula does not hold for class 2!"
    assert abs(qlm[0,2]-stm[0,2]*lambda3)<1e-12, "MMAPPH1PRPR: Little formula does not hold for class 3!"
    
    # ============================= MMAPml.matrix([[K]])/PHml.matrix([[K]])/1 NP tests ============================
    
    print("----------------------------------------------------------------------------")
    
    print("Input:")
    print("------")

    D0=ml.matrix([[-5.49, 0, 1.15, 0], [ 0, -2.29, 0, 0], [ 0, 0.08, -1.32, 0], [ 0.72, 1.17, 0.7, -7.07]])
    D1=ml.matrix([[0.25, 0.38, 0.64, 0], [ 0, 0, 0, 1.09], [ 0, 1.24, 0, 0], [ 0.37, 0, 0, 0]])
    D2=ml.matrix([[0.3, 1.0, 0, 0.48], [ 0, 0.2, 0, 0], [ 0, 0, 0, 0], [ 0.61, 0, 0, 0.2]])
    D3=ml.matrix([[0, 0.98, 0, 0.31], [ 0, 0, 1.0, 0], [ 0, 0, 0, 0], [ 1.1, 0.84, 0.33, 1.03]])
    
    sigma3 = ml.matrix([[1/6, 5/6]])
    S3 = ml.matrix([[-0.5, 0.5], [ 0, -3]])
    sigma2 = ml.matrix([[1/1.7, 1-1/1.7]])
    S2 = ml.matrix([[-2/0.85, 2/0.85], [ 0, -4]])
    sigma1 = ml.matrix([[1/4, 1-1/4]])
    S1 = ml.matrix([[-2.5, 2.5], [ 0, -10]])
    
    print("Test:")
    print("-----")

    qlm, qld = MMAPPH1NPPR((D0,D1,D2,D3), (sigma1,sigma2,sigma3),(S1,S2,S3), 'qlMoms', 3, 'qlDistr', 500)
    momFromDistr = ml.vstack((ml.matrix(np.arange(500))*qld,ml.matrix(np.arange(500)**2)*qld,ml.matrix(np.arange(500)**3)*qld))
    assert la.norm((momFromDistr-qlm)/qlm)<0.001, "MMAPPH1NPPR: queue length moments and queue length distribution are not consistent!"
    
    stm, std = MMAPPH1NPPR((D0,D1,D2,D3), (sigma1,sigma2,sigma3),(S1,S2,S3), 'stMoms', 3, 'stDistr', [1,5,10])   
    assert np.min(std)>=0 and np.max(std)<=1 and np.all(np.diff(std)>=0), "MMAPPH1NPPR: invalid sojourn time distribution!"

    # check Little formula for every traffic class
    lambda1 = 1/MarginalMomentsFromMAP(D0+D2+D3,D1,1)[0]
    lambda2 = 1/MarginalMomentsFromMAP(D0+D1+D3,D2,1)[0]
    lambda3 = 1/MarginalMomentsFromMAP(D0+D1+D2,D3,1)[0]
    assert abs(qlm[0,0]-stm[0,0]*lambda1)<1e-12, "MMAPPH1NPPR: Little formula does not hold for class 1!"
    assert abs(qlm[0,1]-stm[0,1]*lambda2)<1e-12, "MMAPPH1NPPR: Little formula does not hold for class 2!"
    assert abs(qlm[0,2]-stm[0,2]*lambda3)<1e-12, "MMAPPH1NPPR: Little formula does not hold for class 3!"
       
    # ============================= Fluid tests ===============================
    
    print("----------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    Q = ml.matrix([[-9, 2, 4, 0, 1, 2], [ 6, -25, 5, 3, 7, 4], [ 1, 3, -4, 0, 0, 0], [ 0, 0, 0, -8, 3, 5], [ 7, 3, 0, 2, -13, 1], [ 7, 8, 0, 3, 8, -26]])
    Rin = Diag(ml.matrix([[4, 2, 1, 0, 0, 3]]))
    Rout = Diag(ml.matrix([[6, 2, 0, 0, 3, 2]]))
    lambd = np.sum(CTMCSolve(Q)*Rin)
   
    print("Test:")
    print("-----")

    print("qld,qlm = FluidQueue(Q, Rin, Rout, 'qlDistr', np.arange(0,1.1,0.1), 'qlMoms', 5):")
    qld,qlm = FluidQueue(Q, Rin, Rout, 'qlDistr', np.arange(0,1.1,0.1), 'qlMoms', 5)
    print("qld=",qld,"\nqlm=",qlm)
    print("std, stm = FluidQueue(Q, Rin, Rout, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5):")
    std, stm = FluidQueue(Q, Rin, Rout, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5)
    print("std=",std,"\nstm=",stm)
    print("alphap,Ap = FluidQueue(Q, Rin, Rout, 'qlDistrPH'):")
    alphap,Ap = FluidQueue(Q, Rin, Rout, 'qlDistrPH')
    print("alphap=",alphap,",Ap=",Ap)
    print("alpha,A = FluidQueue(Q, Rin, Rout, 'qlDistrME'):")
    alpha,A = FluidQueue(Q, Rin, Rout, 'qlDistrME')
    print("alpha=",alpha,",A=",A)
    print("betap,Bp = FluidQueue(Q, Rin, Rout, 'stDistrPH'):")
    betap, Bp = FluidQueue(Q, Rin, Rout, 'stDistrPH')
    print("betap=",betap,",Bp=",Bp)
    print("beta,B = FluidQueue(Q, Rin, Rout, 'stDistrME'):")
    beta, B = FluidQueue(Q, Rin, Rout, 'stDistrME')
    print("beta=",beta,",B=",B)

    assert CheckMERepresentation(alpha,A), "FluidQueue: invalid ME representation of the queue length!"
    assert CheckMERepresentation(beta,B), "FluidQueue: invalid ME representation of the sojourn time!"
    assert CheckPHRepresentation(alphap,Ap), "FluidQueue: invalid PH representation of the queue length!"
    assert CheckPHRepresentation(betap,Bp), "FluidQueue: invalid PH representation of the sojourn time!"

    # check Little formula
    mql = MomentsFromME(alpha,A,1)[0]
    mst = MomentsFromME(beta,B,1)[0]
    assert abs(mql-mst*lambd)<1e-12, "FluidQueue: Little formula does not hold!"

    # check the equality of the PH and ME results
    assert la.norm((np.array(MomentsFromPH(alphap,Ap,5))-np.array(MomentsFromME(alpha,A,5)))/np.array(MomentsFromME(alpha,A,5)))<1e-12, "FluidQueue: the ME and PH representations are not equal!"
    assert la.norm((np.array(MomentsFromPH(betap,Bp,5))-np.array(MomentsFromME(beta,B,5)))/np.array(MomentsFromME(beta,B,5)))<1e-12, "FluidQueue: the ME and PH representations are not equal!"
    
    # check moment and distribution calculation
    assert la.norm(qld-CdfFromME(alpha,A,np.arange(0,1.1,0.1)))<1e-12, "FluidQueue: qlDistr returns wrong queue length distribution!"
    assert la.norm(std-CdfFromME(beta,B,np.arange(0,1.1,0.1)))<1e-12, "FluidQueue: stDistr returns wrong sojourn time distribution!"
    assert la.norm(np.array(qlm)-np.array(MomentsFromME(alpha,A,5)))<1e-7, "FluidQueue: qlMoms returns wrong queue length moments!"
    assert la.norm(np.array(stm)-np.array(MomentsFromME(beta,B,5)))<1e-7, "FluidQueue: stMoms returns wrong sojourn time moments!"

    # ============================= FluFlu tests ===============================
    
    print("----------------------------------------------------------------------------")
    
    print("Input:")
    print("------")
    
    Qin = ml.matrix([[-2, 1, 1], [ 2, -5, 3], [ 4, 0, -4]])
    Rin = Diag(ml.matrix([[3, 7, 0]]))
    lambd = np.sum(CTMCSolve(Qin)*Rin)
    
    Qout = ml.matrix([[-4, 1, 3], [ 6, -8, 2], [ 3, 7, -10]])
    Rout = Diag(ml.matrix([[1, 7, 15]]))

    print("Test:")
    print("-----")

    print("qld,qlm = FluFluQueue(Qin, Rin, Qout, Rout, False, 'qlDistr', np.arange(0,1.1,0.1), 'qlMoms', 5):")
    qld,qlm = FluFluQueue(Qin, Rin, Qout, Rout, False, 'qlDistr', np.arange(0,1.1,0.1), 'qlMoms', 5)
    print("qld=",qld,"\nqlm=",qlm)
    print("std, stm = FluFluQueue(Qin, Rin, Qout, Rout, False, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5):")
    std, stm = FluFluQueue(Qin, Rin, Qout, Rout, False, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5)
    print("std=",std,"\nstm=",stm)
    print("alphap,Ap = FluFluQueue(Qin, Rin, Qout, Rout, False, 'qlDistrPH'):")
    alphap,Ap = FluFluQueue(Qin, Rin, Qout, Rout, False, 'qlDistrPH')
    print("alphap=",alphap,",Ap=",Ap)
    print("alpha,A = FluFluQueue(Qin, Rin, Qout, Rout, False, 'qlDistrME'):")
    alpha,A = FluFluQueue(Qin, Rin, Qout, Rout, False, 'qlDistrME')
    print("alpha=",alpha,",A=",A)
    print("betap,Bp = FluFluQueue(Qin, Rin, Qout, Rout, False, 'stDistrPH'):")
    betap, Bp = FluFluQueue(Qin, Rin, Qout, Rout, False, 'stDistrPH')
    print("betap=",betap,",Bp=",Bp)
    print("beta,B = FluFluQueue(Qin, Rin, Qout, Rout, False, 'stDistrME'):")
    beta, B = FluFluQueue(Qin, Rin, Qout, Rout, False, 'stDistrME')
    print("beta=",beta,",B=",B)

    assert CheckMERepresentation(alpha,A), "FluFluQueue: invalid ME representation of the queue length!"
    assert CheckMERepresentation(beta,B), "FluFluQueue: invalid ME representation of the sojourn time!"
    assert CheckPHRepresentation(alphap,Ap), "FluFluQueue: invalid PH representation of the queue length!"
    assert CheckPHRepresentation(betap,Bp), "FluFluQueue: invalid PH representation of the sojourn time!"

    # cross-check
    Iin = ml.eye(Qin.shape[0])
    Iout = ml.eye(Qout.shape[0])
    gamma, G = FluidQueue (np.kron(Qin,Iout)+np.kron(Iin,Qout), np.kron(Rin,Iout), np.kron(Iin,Rout), "stDistrME")
    msmall = np.array(MomentsFromME(beta,B,5))
    mlarge = np.array(MomentsFromME(gamma,G,5))
    assert la.norm((msmall-mlarge)/msmall)<1e-12, "FluFluQueue: Large and small model does not give the same results!"

    # check Little formula
    mql = MomentsFromME(alpha,A,1)[0]
    mst = MomentsFromME(beta,B,1)[0]
    assert abs(mql-mst*lambd)<1e-12, "FluFluQueue: Little formula does not hold!"

    # check the equality of the PH and ME results
    assert la.norm((np.array(MomentsFromPH(alphap,Ap,5))-np.array(MomentsFromME(alpha,A,5)))/np.array(MomentsFromME(alpha,A,5)))<1e-12, "FluFluQueue: the ME and PH representations are not equal!"
    assert la.norm((np.array(MomentsFromPH(betap,Bp,5))-np.array(MomentsFromME(beta,B,5)))/np.array(MomentsFromME(beta,B,5)))<1e-12, "FluFluQueue: the ME and PH representations are not equal!"
    
    # check moment and distribution calculation
    assert la.norm(qld-CdfFromME(alpha,A,np.arange(0,1.1,0.1)))<1e-12, "FluFluQueue: qlDistr returns wrong queue length distribution!"
    assert la.norm(std-CdfFromME(beta,B,np.arange(0,1.1,0.1)))<1e-12, "FluFluQueue: stDistr returns wrong sojourn time distribution!"
    assert la.norm(np.array(qlm)-np.array(MomentsFromME(alpha,A,5)))<1e-7, "FluFluQueue: qlMoms returns wrong queue length moments!"
    assert la.norm(np.array(stm)-np.array(MomentsFromME(beta,B,5)))<1e-7, "FluFluQueue: stMoms returns wrong sojourn time moments!"


    print("Test:")
    print("-----")

    print("qld,qlm = FluFluQueue(Qin, Rin, Qout, Rout, True, 'qlDistr', np.arange(0,1.1,0.1), 'qlMoms', 5):")
    qld,qlm = FluFluQueue(Qin, Rin, Qout, Rout, True, 'qlDistr', np.arange(0,1.1,0.1), 'qlMoms', 5)
    print("qld=",qld,"\nqlm=",qlm)
    print("std, stm = FluFluQueue(Qin, Rin, Qout, Rout, True, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5):")
    std, stm = FluFluQueue(Qin, Rin, Qout, Rout, True, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5)
    print("std=",std,"\nstm=",stm)
    print("alphap,Ap = FluFluQueue(Qin, Rin, Qout, Rout, True, 'qlDistrPH'):")
    alphap,Ap = FluFluQueue(Qin, Rin, Qout, Rout, True, 'qlDistrPH')
    print("alphap=",alphap,",Ap=",Ap)
    print("alpha,A = FluFluQueue(Qin, Rin, Qout, Rout, True, 'qlDistrME'):")
    alpha,A = FluFluQueue(Qin, Rin, Qout, Rout, True, 'qlDistrME')
    print("alpha=",alpha,",A=",A)
    print("betap,Bp = FluFluQueue(Qin, Rin, Qout, Rout, True, 'stDistrPH'):")
    betap, Bp = FluFluQueue(Qin, Rin, Qout, Rout, True, 'stDistrPH')
    print("betap=",betap,",Bp=",Bp)
    print("beta,B = FluFluQueue(Qin, Rin, Qout, Rout, True, 'stDistrME'):")
    beta, B = FluFluQueue(Qin, Rin, Qout, Rout, True, 'stDistrME')
    print("beta=",beta,",B=",B)

    assert CheckMERepresentation(alpha,A), "FluFluQueue: invalid ME representation of the queue length!"
    assert CheckMERepresentation(beta,B), "FluFluQueue: invalid ME representation of the sojourn time!"
    assert CheckPHRepresentation(alphap,Ap), "FluFluQueue: invalid PH representation of the queue length!"
    assert CheckPHRepresentation(betap,Bp), "FluFluQueue: invalid PH representation of the sojourn time!"

    # cross-check
    Iin = ml.eye(Qin.shape[0])
    Iout = ml.eye(Qout.shape[0])
    gamma, G = FluidQueue (np.kron(Qin,Iout)+np.kron(Iin,Qout), np.kron(Rin,Iout), np.kron(Iin,Rout), 'Q0', np.kron(Qin,Iout)+np.kron(Rin,la.pinv(Rout)*Qout), "stDistrME")
    msmall = np.array(MomentsFromME(beta,B,5))
    mlarge = np.array(MomentsFromME(gamma,G,5))
    assert la.norm((msmall-mlarge)/msmall)<1e-12, "FluFluQueue: Large and small model does not give the same results!"

    # check Little formula
    mql = MomentsFromME(alpha,A,1)[0]
    mst = MomentsFromME(beta,B,1)[0]
    assert abs(mql-mst*lambd)<1e-12, "FluFluQueue: Little formula does not hold!"

    # check the equality of the PH and ME results
    assert la.norm((np.array(MomentsFromPH(alphap,Ap,5))-np.array(MomentsFromME(alpha,A,5)))/np.array(MomentsFromME(alpha,A,5)))<1e-12, "FluFluQueue: the ME and PH representations are not equal!"
    assert la.norm((np.array(MomentsFromPH(betap,Bp,5))-np.array(MomentsFromME(beta,B,5)))/np.array(MomentsFromME(beta,B,5)))<1e-12, "FluFluQueue: the ME and PH representations are not equal!"
    
    # check moment and distribution calculation
    assert la.norm(qld-CdfFromME(alpha,A,np.arange(0,1.1,0.1)))<1e-12, "FluFluQueue: qlDistr returns wrong queue length distribution!"
    assert la.norm(std-CdfFromME(beta,B,np.arange(0,1.1,0.1)))<1e-12, "FluFluQueue: stDistr returns wrong sojourn time distribution!"
    assert la.norm(np.array(qlm)-np.array(MomentsFromME(alpha,A,5)))<1e-7, "FluFluQueue: qlMoms returns wrong queue length moments!"
    assert la.norm(np.array(stm)-np.array(MomentsFromME(beta,B,5)))<1e-7, "FluFluQueue: stMoms returns wrong sojourn time moments!"

if __name__ == "__main__":
    TestQueuesPackage()