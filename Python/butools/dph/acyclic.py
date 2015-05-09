import numpy as np
import butools
from butools.reptrans import SimilarityMatrix
from butools.dph import CheckDPHRepresentation, CheckMGRepresentation
import scipy.linalg as la
import numpy.matlib as ml

def AcyclicDPHFromMG (alpha, A, maxSize=100, precision=1e-14):
    """
    Transforms a matrix-geometric representation to an acyclic
    DPH representation of the same size, if possible.
    
    Parameters
    ----------
    alpha : matrix, shape (1,N)
        Initial vector of the distribution
    A : matrix, shape (N,N)
        Matrix parameter of the distribution
    precision : double, optional
        Vector and matrix entries smaller than the precision
        are considered to be zeros. The default value is 1e-14.
    
    Returns
    -------
    beta : matrix, shape (1,M)
        The initial probability vector of the acyclic discrete
        phase-type representation
    B : matrix, shape (M,M)
        Transition probability matrix of the acyclic discrete
        phase-type representation
    
    Notes
    -----
    Contrary to 'AcyclicPHFromME' of the 'ph' package, this 
    procedure is not able to extend the size in order to obtain
    a Markovian initial vector.
    
    Raises an error if A has complex eigenvalues. In this case
    the transformation to an acyclic representation is not 
    possible
    """

    if butools.checkInput and not CheckMGRepresentation(alpha,A):
        raise Exception("AcyclicDPHFromMG: input is not a valid MG representation!")

    ev = la.eigvals(A)
    ix = np.argsort(np.abs(np.real(ev)))
    lambda2 = ev[ix]
    lambda3 = lambda2[lambda2!=lambda2[-1]]

    if np.max(np.abs(np.imag(ev)))>precision:
        raise Exception("AcyclicDPHFromMG: The input matrix has complex eigenvalue!")

    mx = ml.matrix(np.diag(lambda2)+np.diag(1.0-lambda3, 1))
    T = SimilarityMatrix (A, mx)
    gamma = alpha*T
    G = mx
    if not CheckDPHRepresentation (gamma, G, precision):
        raise Exception("AcyclicDPHFromMG: No acyclic representation found!")
    else:
        return (gamma, G)
