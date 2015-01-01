from numpy import linalg as la
import numpy.matlib as ml
import numpy as np
from butools.moments import ReducedMomsFromMoms, FactorialMomsFromMoms
from butools.ph import MEFromMoments

def MGFromMoments (moms):

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
