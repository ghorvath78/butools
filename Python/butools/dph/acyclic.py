import numpy as np
import butools
from butools.reptrans import SimilarityMatrix
from butools.dph import CheckDPHRepresentation, CheckMGRepresentation
import scipy.linalg as la
import numpy.matlib as ml

def AcyclicDPHFromMG (alpha, A, maxSize=100, precision=1e-14):

    if butools.checkInput and not CheckMGRepresentation(alpha,A,precision):
        raise Exception("AcyclicDPHFromMG: input is not a valid MG representation!")

    ev = la.eigvals(A)
    ix = np.argsort(np.abs(np.real(ev)))
    lambda2 = ev[ix]
    lambda3 = lambda2[lambda2!=lambda2[-1]]

    mx = ml.matrix(np.diag(lambda2)+np.diag(1.0-lambda3, 1))
    T = SimilarityMatrix (A, mx)
    gamma = alpha*T
    G = mx
    if not CheckDPHRepresentation (gamma, G, precision):
        raise Exception("AcyclicDPHFromMG: No acyclic representation found!")
    else:
        return (gamma, G)
