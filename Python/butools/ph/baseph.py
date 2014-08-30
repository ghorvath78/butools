import numpy as np
from numpy import linalg as la
from scipy.linalg import expm
import numpy.matlib as ml
import math
import butools
from butools.ph import CheckMERepresentation, CheckPHRepresentation
from butools.mc import CTMCSolve
from butools.utils import Diag

def MomentsFromME (alpha, A, K=0, prec=1e-14):

    if butools.checkInput and not CheckMERepresentation (alpha, A, prec):
        raise Exception("MomentsFromME: Input is not a valid ME representation!")
        
    if K==0:
        K = 2*np.size(alpha,1)-1
    return [math.factorial(i)*np.sum(alpha*la.inv(-A)**i) for i in range(1,K+1)]
  
def MomentsFromPH (alpha, A, K=0, prec=1e-14):

    if butools.checkInput and not CheckPHRepresentation (alpha, A, prec):
        raise Exception("MomentsFromPH: Input is not a valid PH representation!")

    return MomentsFromME (alpha, A, K, prec)

def PdfFromME (alpha, A, x, prec=1e-14):

    if butools.checkInput and not CheckMERepresentation (alpha, A, prec):
        raise Exception("PdfFromME: Input is not a valid ME representation!")

    y = [np.sum(alpha*expm(A*xv)*(-A)) for xv in x]
    return np.array(y)

def PdfFromPH (alpha, A, x, prec=1e-14):

    if butools.checkInput and not CheckPHRepresentation (alpha, A, prec):
        raise Exception("PdfFromPH: Input is not a valid PH representation!")

    return PdfFromME (alpha, A, x, prec)

def IntervalPdfFromPH (alpha, A, intBounds, prec=1e-14):

    if butools.checkInput and not CheckPHRepresentation (alpha, A, prec):
        raise Exception("IntervalPdfFromPH: Input is not a valid PH representation!")

    steps = len(intBounds)
    x = [(intBounds[i+1]+intBounds[i])/2.0 for i in range(steps-1)];
    y = [(np.sum(alpha*expm(A*intBounds[i])) - np.sum(alpha*expm(A*intBounds[i+1])))/(intBounds[i+1]-intBounds[i]) for i in range(steps-1)]
    return (np.array(x), np.array(y))

def CdfFromME (alpha, A, x, prec=1e-14):

    if butools.checkInput and not CheckMERepresentation (alpha, A, prec):
        raise Exception("CdfFromME: Input is not a valid ME representation!")

    return np.array([1.0-np.sum(alpha*expm(A*xv)) for xv in x])

def CdfFromPH (alpha, A, x, prec=1e-14):

    if butools.checkInput and not CheckPHRepresentation (alpha, A, prec):
        raise Exception("CdfFromPH: Input is not a valid PH representation!")

    return CdfFromME (alpha, A, x, prec)
    
def RandomPH (order, mean=1.0, zeroEntries=0, maxTrials=1000, prec=1e-14):

    if zeroEntries > (order+1)*(order-1):
        raise Exception ("RandomPH: Too many zero entries requested! Try to decrease the zero entries number!")

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
            # construct PH parameters
            A = ml.matrix(B[:,:order])
            a = ml.matrix(B[:,order+1]).T
            A = A - Diag(np.sum(A,1)+a)
            alpha = ml.matrix(B[:,order])
            # check if it is a proper PH (irreducible phase process & no full zero matrix)
            if np.all(A==0.0) or np.all(alpha==0.0) or np.all(a==0.0):
                continue
            alpha = alpha / np.sum(alpha)
            D = A + a*alpha
            if la.matrix_rank(D) == order-1:
                pi = CTMCSolve(D, prec)
                if np.min(np.abs(pi)) > math.sqrt(prec):
                    # scale to the mean value
                    m = MomentsFromPH (alpha, A, 1, prec)[0]
                    A *= m / mean
                    return (alpha, A)
            trials += 1
    raise Exception("No feasible random PH found with such many zero entries!")    
    