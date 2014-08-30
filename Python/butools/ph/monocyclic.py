import numpy as np
from butools.reptrans import TransformToMonocyclic, TransformToAcyclic, SimilarityMatrix, ExtendToMarkovian
from butools.ph import CheckPHRepresentation, CheckMERepresentation

def MonocyclicPHFromME (alpha, A, maxSize=100, precision=1e-14):
        
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

def CheckMEPositiveDensity (alpha, A, maxSize=100, prec=1e-14):

    try:
        beta, B = MonocyclicPHFromME (alpha, A, maxSize, prec)
        r = CheckMERepresentation (beta, B, prec)
    except Exception:
        r = False
    return r
    