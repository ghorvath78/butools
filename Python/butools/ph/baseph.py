import numpy as np
from numpy import linalg as la
from scipy.linalg import expm, expm2
import numpy.matlib as ml
import math
import butools
from butools.ph import CheckMERepresentation, CheckPHRepresentation
from butools.mc import CTMCSolve
from butools.utils import Diag

def MomentsFromME (alpha, A, K=0):
    """
    Returns the first K moments of a matrix-exponential
    distribution.
    
    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-exponential
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential
        distribution.
    K : int, optional
        Number of moments to compute. If K=0, 2*M-1 moments
        are computed. The default value is K=0.
    prec : double, optional
        Numerical precision for checking the input.
        The default value is 1e-14.
    
    Returns
    -------
    moms : row vector of doubles
        The vector of moments.
        
    """

    if butools.checkInput and not CheckMERepresentation (alpha, A):
        raise Exception("MomentsFromME: Input is not a valid ME representation!")
        
    if K==0:
        K = 2*np.size(alpha,1)-1
    return [math.factorial(i)*np.sum(alpha*la.inv(-A)**i) for i in range(1,K+1)]
  
def MomentsFromPH (alpha, A, K=0):
    """
    Returns the first K moments of a continuous phase-type
    distribution.
    
    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial probability vector of the phase-type
        distribution.
    A : matrix, shape (M,M)
        The transient generator matrix of the phase-type
        distribution.
    K : int, optional
        Number of moments to compute. If K=0, 2*M-1 moments
        are computed. The default value is K=0.
    prec : double, optional
        Numerical precision for checking the input.
        The default value is 1e-14.
    
    Returns
    -------
    moms : row vector of doubles
        The vector of moments.
    """

    if butools.checkInput and not CheckPHRepresentation (alpha, A):
        raise Exception("MomentsFromPH: Input is not a valid PH representation!")

    return MomentsFromME (alpha, A, K)

def PdfFromME (alpha, A, x):
    """
    Returns the probability density function of a matrix-
    exponential distribution.
    
    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-exponential
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential
        distribution.
    x : vector of doubles
        The density function will be computed at these points
    prec : double, optional
        Numerical precision to check if the input ME 
        distribution is valid. The default value is 1e-14.
    
    Returns
    -------
    pdf : column vector of doubles
        The values of the density function at the 
        corresponding "x" values
    """

    if butools.checkInput and not CheckMERepresentation (alpha, A):
        raise Exception("PdfFromME: Input is not a valid ME representation!")

    y = [np.sum(alpha*expm(A*xv)*(-A)) for xv in x]
    return np.array(y)

def PdfFromPH (alpha, A, x):
    """
    Returns the probability density function of a continuous
    phase-type distribution.
    
    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial probability vector of the phase-type
        distribution.
    A : matrix, shape (M,M)
        The transient generator matrix of the phase-type
        distribution.
    x : vector of doubles
        The density function will be computed at these points
    prec : double, optional
        Numerical precision to check if the input ME 
        distribution is valid. The default value is 1e-14.
    
    Returns
    -------
    pdf : column vector of doubles
        The values of the density function at the 
        corresponding "x" values
    """

    if butools.checkInput and not CheckPHRepresentation (alpha, A):
        raise Exception("PdfFromPH: Input is not a valid PH representation!")

    return PdfFromME (alpha, A, x)

def IntervalPdfFromPH (alpha, A, intBounds):
    """
    Returns the approximate probability density function of a
    continuous phase-type distribution, based on the 
    probability of falling into intervals.
    
    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial probability vector of the phase-type
        distribution.
    A : matrix, shape (M,M)
        The transient generator matrix of the phase-type
        distribution.
    intBounds : vector, shape (K)
        The array of interval boundaries. The pdf is the
        probability of falling into an interval divided by
        the interval length. 
        If the size of intBounds is K, the size of the result is K-1.
    prec : double, optional
        Numerical precision to check if the input is a valid
        phase-type distribution. The default value is 1e-14
    
    Returns
    -------
    x : matrix of doubles, shape(K-1,1)
        The points at which the pdf is computed. It holds the center of the 
        intervals defined by intBounds.
    y : matrix of doubles, shape(K-1,1)
        The values of the density function at the corresponding "x" values
    
    Notes
    -----
    This method is more suitable for comparisons with empirical
    density functions than the exact one (given by PdfFromPH).
    """

    if butools.checkInput and not CheckPHRepresentation (alpha, A):
        raise Exception("IntervalPdfFromPH: Input is not a valid PH representation!")

    steps = len(intBounds)
    x = [(intBounds[i+1]+intBounds[i])/2.0 for i in range(steps-1)]
    y = [(np.sum(alpha*expm2((A*intBounds[i]).A)) - np.sum(alpha*expm2((A*intBounds[i+1]).A)))/(intBounds[i+1]-intBounds[i]) for i in range(steps-1)]
    return (np.array(x), np.array(y))

def CdfFromME (alpha, A, x):
    """
    Returns the cummulative distribution function of a
    matrix-exponential distribution.
    
    Parameters
    ----------
    alpha : matrix, shape (1,M)
        The initial vector of the matrix-exponential
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential
        distribution.
    x : vector of doubles
        The cdf will be computed at these points
    
    Returns
    -------
    cdf : column vector of doubles
        The values of the cdf at the corresponding "x" values
    """

    if butools.checkInput and not CheckMERepresentation (alpha, A):
        raise Exception("CdfFromME: Input is not a valid ME representation!")

    return np.array([1.0-np.sum(alpha*expm(A.A*xv)) for xv in x])

def CdfFromPH (alpha, A, x):
    """
    Returns the cummulative distribution function of a
    continuous phase-type distribution.
    
    Parameters
    ----------
    alpha : matrix, shape (1,M)
        The initial probability vector of the phase-type
        distribution.
    A : matrix, shape (M,M)
        The transient generator matrix of the phase-type
        distribution.
    x : vector of doubles
        The cdf will be computed at these points
    
    Returns
    -------
    cdf : column vector of doubles
        The values of the cdf at the corresponding "x" values
    """

    if butools.checkInput and not CheckPHRepresentation (alpha, A):
        raise Exception("CdfFromPH: Input is not a valid PH representation!")

    return CdfFromME (alpha, A, x)
    
def RandomPH (order, mean=1.0, zeroEntries=0, maxTrials=1000, prec=1e-7):
    """
    Returns a random phase-type distribution with a given 
    order.
    
    Parameters
    ----------
    order : int
        The size of the phase-type distribution
    mean : double, optional
        The mean of the phase-type distribution 
    zeroEntries : int, optional
        The number of zero entries in the initial vector, 
        generator matrix and closing vector
    maxTrials : int, optional
        The maximum number of trials to find a proper PH 
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
    parameter.   
    """

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
                pi = CTMCSolve(D)
                if np.min(np.abs(pi)) > prec:
                    # scale to the mean value
                    m = MomentsFromPH (alpha, A, 1)[0]
                    A *= m / mean
                    return (alpha, A)
            trials += 1
    raise Exception("No feasible random PH found with such many zero entries!")    
    
