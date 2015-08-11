import numpy as np
from numpy import linalg as la
import numpy.matlib as ml
import math
import butools
from butools.dph import CheckMGRepresentation, CheckDPHRepresentation
from butools.moments import MomsFromFactorialMoms
from butools.mc import DTMCSolve
from butools.utils import Diag


def MomentsFromMG (alpha, A, K=0):
    """
    Returns the first K moments of a matrix geometric 
    distribution.
    
    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-geometric distribution.
        The sum of the entries of alpha is less or equal to 1.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-geometric 
        distribution.
    K : int, optional
        Number of moments to compute. If K=0, 2*M-1 moments are
        computed. The default value is 0.
    prec : double, optional
        Numerical precision for checking the input.
        The default value is 1e-14.
    
    Returns
    -------
    moms : row vector of doubles
        The vector of moments.
        
    """

    if butools.checkInput and not CheckMGRepresentation (alpha, A):
        raise Exception("MomentsFromMG: Input is not a valid MG representation!")

    m = A.shape[0]
    if K==0:
        K = 2*m-1
    Ai = la.inv(ml.eye(m)-A)
    return MomsFromFactorialMoms([math.factorial(i)*np.sum(alpha*Ai**i*A**(i-1)) for i in range(1,K+1)])

def MomentsFromDPH (alpha, A, K=0):
    """
    Returns the first K moments of a discrete phase-type
    distribution.
    
    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial probability vector of the discrete phase-
        type distribution. The sum of the entries of pi0 is 
        less or equal to 1.
    A : matrix, shape (M,M)
        The transient generator matrix of the discrete phase-
        type distribution.
    K : int, optional
        Number of moments to compute. If K=0, 2*M-1 moments
        are computed. The default value is 0.
    prec : double, optional
        Numerical precision for checking the input.
        The default value is 1e-14.
    
    Returns
    -------
    moms : row vector of doubles
        The vector of moments.
        
    """

    if butools.checkInput and not CheckDPHRepresentation (alpha, A):
        raise Exception("MomentsFromDPH: Input is not a valid DPH representation!")

    return MomentsFromMG (alpha, A, K)

def PmfFromMG (alpha, A, x):
    """
    Returns the probability mass function of a matrix-
    geometric distribution.
    
    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-geometric
        distribution. The sum of the entries of pi0 is less
        or equal to 1.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-geometric
        distribution.
    x : vector of non-negative integers
        The density function will be computed at these points
    prec : double, optional
        Numerical precision to check if the input MG 
        distribution is valid. The default value is 1e-14.
    
    Returns
    -------
    pmf : column vector of doubles
        The probabilities that the matrix-geometrically 
        distributed random variable takes the corresponding "x"
        values
        
    """

    if butools.checkInput and not CheckMGRepresentation (alpha, A):
        raise Exception("PmfFromMG: Input is not a valid MG representation!")

    a = 1-np.sum(A,1)
    y = np.empty(len(x))
    for i in range(len(y)):
        if x[i]==0:
            y[i] = 1.0 - np.sum(alpha)
        else:
            y[i] = np.sum(alpha*(A**int(x[i]-1))*a)
    return y

def PmfFromDPH (alpha, A, x):
    """
    Returns the probability mass function of a discrete
    phase-type distribution.
    
    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial probability vector of the discrete phase-
        type distribution. The sum of the entries of pi0 is
        less or equal to 1.
    A : matrix, shape (M,M)
        The transient generator matrix of the discrete phase-
        type distribution.
    x : vector of non-negative integers
        The density function will be computed at these points
    prec : double, optional
        Numerical precision to check if the input DPH 
        distribution is valid. The default value is 1e-14.
    
    Returns
    -------
    pmf : column vector of doubles
        The probabilities that the discrete phase type
        distributed random variable takes the corresponding
        "x" values
        
    """

    if butools.checkInput and not CheckDPHRepresentation (alpha, A):
        raise Exception("PmfFromDPH: Input is not a valid DPH representation!")

    return PmfFromMG (alpha, A, x)

def CdfFromMG (alpha, A, x):
    """
    Returns the cummulative distribution function of a 
    matrix-geometric distribution.
    
    Parameters
    ----------
    alpha : matrix, shape (1,M)
        The initial vector of the matrix-geometric distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-geometric 
        distribution.
    x : vector of non-negative integers
        The density function will be computed at these points
    prec : double, optional
        Numerical precision to check if the input MG 
        distribution is valid. The default value is 1e-14.
    
    Returns
    -------
    cdf : column vector of doubles
        The probabilities that the matrix-geometrically 
        distributed random variable is less or equal to
        the corresponding "x" values
        
    """

    if butools.checkInput and not CheckMGRepresentation (alpha, A):
        raise Exception("CdfFromMG: Input is not a valid MG representation!")

    y = np.empty(len(x))
    for i in range(len(y)):
        y[i] = 1.0-np.sum(alpha*(A**int(x[i])))
    return y

def CdfFromDPH (alpha, A, x):
    """
    Returns the cummulative distribution function of a 
    discrete phase-type distribution.
    
    Parameters
    ----------
    alpha : matrix, shape (1,M)
        The initial probability vector of the discrete phase-
        type distribution.
    A : matrix, shape (M,M)
        The transition probability  matrix of the discrete phase-
        type distribution.
    x : vector of non-negative integers
        The cdf will be computed at these points
    prec : double, optional
        Numerical precision to check if the input DPH 
        distribution is valid. The default value is 1e-14.
    
    Returns
    -------
    cdf : column vector of doubles
        The probabilities that the discrete phase type 
        distributed random variable is less or equal to the
        corresponding "x" values
        
    """

    if butools.checkInput and not CheckDPHRepresentation (alpha, A):
        raise Exception("CdfFromDPH: Input is not a valid DPH representation!")

    return CdfFromMG (alpha, A, x)

def RandomDPH (order, mean=10.0, zeroEntries=0, maxTrials=1000, prec=1e-7):
    """
    Returns a random discrete phase-type distribution with a 
    given mean value.
    
    Parameters
    ----------
    order : int
        The size of the discrete phase-type distribution
    mean : double, optional
        The mean of the discrete phase-type distribution 
    zeroEntries : int, optional
        The number of zero entries in the initial vector, 
        generator matrix and closing vector
    maxTrials : int, optional
        The maximum number of trials to find a proper DPH 
        (that has an irreducible phase process and none of 
        its parameters is all-zero). The default value is 
        1000.
    prec : double, optional
        Numerical precision for checking the irreducibility.
        The default value is 1e-14.
    
    Returns
    -------
    alpha : vector, shape (1,M)
        The initial probability vector of the phase-type 
        distribution.
    A : matrix, shape (M,M)
        The transient generator matrix of the phase-type 
        distribution.
    
    Notes
    -----
    If the procedure fails, try to increase the 'maxTrials'
    parameter, or increase the mean value.
    """

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
            if la.matrix_rank(ml.eye(order)-A) == order:
                if np.min(np.abs(alpha*la.inv(ml.eye(order)-A))) > prec:
                    # diagonals of matrix A:
                    d = np.random.rand(order)
                    # scale to the mean value
                    m = MomentsFromDPH (alpha, Diag(1-d)*A+Diag(d), 1)[0]
                    d = 1 - (1-d)*m/mean
                    A = Diag(1-d)*A+Diag(d)
                    if CheckDPHRepresentation(alpha,A,prec):
                        return (alpha, A)
            trials += 1
    raise Exception("No feasible random PH found with such many zero entries!")    
    