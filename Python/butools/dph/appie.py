from numpy import linalg as la
import numpy.matlib as ml
import numpy as np
from butools.moments import ReducedMomsFromMoms, FactorialMomsFromMoms
from butools.ph import MEFromMoments

def MGFromMoments (moms):
    """
    Creates a matrix-geometric distribution that has the
    same moments as given.
    
    Parameters
    ----------
    moms : vector of doubles
        The list of moments. The order of the resulting 
        matrix-geometric distribution is 
        determined based on the number of moments given. To 
        obtain a matrix-geometric distribution of order M,
        2*M-1 moments are required.
    
    Returns
    -------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-geometric 
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-geometric 
        distribution.
    
    References
    ----------
    .. [1] A. van de Liefvoort. The moment problem for 
           continuous distributions. Technical report, 
           University of Missouri, WP-CM-1990-02, Kansas City,
           1990.
    """

    rfmoms = ReducedMomsFromMoms (FactorialMomsFromMoms(moms))
    rfmoms[:0] = [1]

    vlist = np.zeros(len(moms))
    tmpVec= np.zeros(len(moms)+1)
    k=1
    tmpVec[0]=rfmoms[0]
    for i in range(len(moms)):
        tmpVec[i+1]=(-1)**(i+1)*rfmoms[i+1]
        k=k*(i+1)
        vlist[i]=k*sum(tmpVec)

    alpha,C = MEFromMoments(vlist)
    iC=la.inv(C)
    A=iC*la.inv(iC+ml.eye(len(iC)))
    return alpha, A
