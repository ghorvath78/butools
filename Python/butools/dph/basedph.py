import numpy as np
from numpy import linalg as la
import numpy.matlib as ml
import math
import butools
from butools.dph import CheckMGRepresentation, CheckDPHRepresentation
from butools.moments import MomsFromFactorialMoms
from butools.mc import DTMCSolve
from butools.utils import Diag


def MomentsFromMG (alpha, A, K=0, prec=1e-14):

    if butools.checkInput and not CheckMGRepresentation (alpha, A, prec):
        raise Exception("MomentsFromMG: Input is not a valid MG representation!")

    m = A.shape[0]
    if K==0:
        K = 2*m-1
    Ai = la.inv(ml.eye(m)-A)
    return MomsFromFactorialMoms([math.factorial(i)*np.sum(alpha*Ai**i*A**(i-1)) for i in range(1,K+1)])

def MomentsFromDPH (alpha, A, K=0, prec=1e-14):

    if butools.checkInput and not CheckDPHRepresentation (alpha, A, prec):
        raise Exception("MomentsFromDPH: Input is not a valid DPH representation!")

    return MomentsFromMG (alpha, A, K)

def PmfFromMG (alpha, A, x, prec=1e-14):

    if butools.checkInput and not CheckMGRepresentation (alpha, A, prec):
        raise Exception("PmfFromMG: Input is not a valid MG representation!")

    a = 1-np.sum(A,1)
    y = np.empty(len(x))
    for i in range(len(y)):
        if x[i]==0:
            y[i] = 1.0 - np.sum(alpha)
        else:
            y[i] = np.sum(alpha*(A**int(x[i]-1))*a)
    return y

def PmfFromDPH (alpha, A, x, prec=1e-14):

    if butools.checkInput and not CheckDPHRepresentation (alpha, A, prec):
        raise Exception("PmfFromDPH: Input is not a valid DPH representation!")

    return PmfFromMG (alpha, A, x)

def CdfFromMG (alpha, A, x, prec=1e-14):

    if butools.checkInput and not CheckMGRepresentation (alpha, A, prec):
        raise Exception("CdfFromMG: Input is not a valid MG representation!")

    y = np.empty(len(x))
    for i in range(len(y)):
        y[i] = 1.0-np.sum(alpha*(A**int(x[i])))
    return y

def CdfFromDPH (alpha, A, x, prec=1e-14):

    if butools.checkInput and not CheckDPHRepresentation (alpha, A, prec):
        raise Exception("CdfFromDPH: Input is not a valid DPH representation!")

    return CdfFromMG (alpha, A, x)

def RandomDPH (order, mean=1.0, zeroEntries=0, maxTrials=1000, prec=1e-14):

    if zeroEntries > (order+1)*(order-1):
        raise Exception ("RandomDPH: Too many zero entries requested! Try to decrease the zero entries number!")

    # distribute the zero entries among the rows
    def allZeroDistr (states, zeros):
        if states==1:
            return [[zeros]]
        else:
            o = [];
            for i in range(zeros+1):
                x = allZeroDistr (states-1, zeros-i)
                for j in range(len(x)):
                    xt = x[j]
                    xt.append(i)
                    xt.sort()
                    # check if we have it already
                    if o.count(xt)==0:
                        o.append(xt)
            return o
    zeroDistr = allZeroDistr(order, zeroEntries)   

    trials = 1
    while trials<maxTrials:
        # select a configuration from zeroDistr: it is a list describing the zero entries in each row
        zdix = np.random.permutation(len(zeroDistr))
        for k in range(len(zeroDistr)):
            zDistr = zeroDistr[zdix[k]];
            
            B = np.zeros((order,order+2))
            for i in range(order):
                rp = np.random.permutation(order+1)
                a = np.zeros(order+1)
                for j in range(order+1 - zDistr[i]):
                    a[rp[j]] = np.random.rand()
                B[i,0:i] = a[0:i]
                B[i,i+1:] = a[i:]
            # construct DPH parameters
            A = ml.matrix(B[:,:order])
            a = ml.matrix(B[:,order+1]).T

            sc = np.sum(A,1)+a.A
            if np.any(sc==0):
                continue
            A = Diag(1/sc)*A
            a = Diag(1/sc)*a            
            alpha = ml.matrix(B[:,order])
            # check if it is a proper PH (irreducible phase process & no full zero matrix)
            if np.all(A==0.0) or np.all(alpha==0.0) or np.all(a==0.0):
                continue
            alpha = alpha / np.sum(alpha)
            D = A + a*alpha
            if la.matrix_rank(D) == order-1:
                pi = DTMCSolve(D, prec)
                if np.min(np.abs(pi)) > math.sqrt(prec):
                    # diagonals of matrix A:
                    d = np.random.rand(order)
                    # scale to the mean value
                    m = MomentsFromDPH (alpha, Diag(1-d)*A+Diag(d), 1, prec)[0]
                    d = 1 - (1-d)*m/mean
                    A = Diag(1-d)*A+Diag(d)
                    if CheckDPHRepresentation(alpha,A,prec):
                        return (alpha, A)
            trials += 1
    raise Exception("No feasible random PH found with such many zero entries!")    
    