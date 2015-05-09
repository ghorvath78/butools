import butools
import numpy as np
from butools.reptrans import TransformToMonocyclic, TransformToAcyclic, SimilarityMatrix, ExtendToMarkovian
from butools.ph import CheckPHRepresentation, CheckMERepresentation

def MonocyclicPHFromME (alpha, A, maxSize=100, precision=1e-14):
    """
    Transforms an arbitrary matrix-exponential representation
    to a Markovian monocyclic representation.
    
    Parameters
    ----------
    alpha : matrix, shape (1,N)
        Initial vector of the distribution
    A : matrix, shape (N,N)
        Matrix parameter of the distribution
    maxSize : int, optional
        The maximum number of phases for the result. The default
        value is 100.
    precision : double, optional
        Vector and matrix entries smaller than the precision
        are considered to be zeros. The default value is 1e-14.
    
    Returns
    -------
    beta : matrix, shape (1,M)
        The initial probability vector of the Markovian 
        monocyclic representation
    B : matrix, shape (M,M)
        Transient generator matrix of the Markovian 
        monocyclic representation
    
    Notes
    -----
    Raises an error if no Markovian monocyclic representation
    has been found.
    
    References
    ----------
    .. [1]  Mocanu, S., Commault, C.: "Sparse representations of
           phase-type distributions," Stoch. Models 15, 759-778 
           (1999)
    """
        
    G = TransformToMonocyclic (A, maxSize, precision)

    # find transformation matrix
    T = SimilarityMatrix (A, G)
    gamma = np.real(alpha*T)

    if np.min(gamma) <= -precision:
        gamma, G = ExtendToMarkovian (gamma, G, maxSize, precision)

    if not CheckPHRepresentation (gamma, G, precision):
        raise Exception("MonocyclicPHFromME: No monocyclic representation found up to the given size and precision!")
    else:
        return (gamma, G)

def AcyclicPHFromME (alpha, A, maxSize=100, precision=1e-14):
    """
    Transforms an arbitrary matrix-exponential representation
    to an acyclic phase-type representation. (see [1]_).
    
    Parameters
    ----------
    alpha : matrix, shape (1,N)
        Initial vector of the distribution
    A : matrix, shape (N,N)
        Matrix parameter of the distribution
    maxSize : int, optional
        The maximum number of phases for the result.
        The default value is 100.
    precision : double, optional
        Vector and matrix entries smaller than the precision
        are considered to be zeros. The default value is 1e-14.
    
    Returns
    -------
    beta : matrix, shape (1,M)
        The initial probability vector of the Markovian 
        acyclic representation
    B : matrix, shape (M,M)
        Transient generator matrix of the Markovian 
        acyclic representation
    
    Notes
    -----
    Raises an error if no Markovian acyclic representation
    has been found.
    
    References
    ----------
    .. [1]  Mocanu, S., Commault, C.: "Sparse representations of
            phase-type distributions," Stoch. Models 15, 759-778 
            (1999)
    """
    
    G = TransformToAcyclic (A, maxSize, precision)

    # find transformation matrix
    T = SimilarityMatrix (A, G)
    gamma = np.real(alpha*T)

    if np.min(gamma) <= -precision:
        gamma, G = ExtendToMarkovian (gamma, G, maxSize, precision)
    
    if not CheckPHRepresentation (gamma, G, precision):
        raise Exception("AcyclicPHFromME: No acyclic representation found up to the given size and precision!")
    else:
        return (gamma, G)

def CheckMEPositiveDensity (alpha, A, maxSize=100, prec=None):
    """
    Checks if the given matrix-exponential distribution has 
    positive density.
    
    Parameters
    ----------
    alpha : matrix, shape (1,M)
        Initial vector of the matrix-exponential distribution 
        to check
    A : matrix, shape (M,M)
        Matrix parameter of the matrix-exponential distribution
        to check
    maxSize : int, optional
        The procedure tries to transform the ME distribution
        to phase-type up to order maxSize. The default value
        is 100.
    prec : double, optional
        Numerical precision. The default value is 1e-14.
    
    Returns
    -------
    r : bool
        True, if the given matrix-exponential distribution has
        a positive density
    
    Notes
    -----
    This procedure calls MonocyclicPHFromME, and can be time 
    consuming. 
    """

    if prec==None:
        prec=butools.checkPrecision

    try:
        beta, B = MonocyclicPHFromME (alpha, A, maxSize, prec)
        r = CheckMERepresentation (beta, B, prec)
    except Exception:
        r = False
    return r
    