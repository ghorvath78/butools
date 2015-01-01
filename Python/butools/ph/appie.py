import numpy as np
from numpy import linalg as la
import numpy.matlib as ml
from butools.moments import ReducedMomsFromMoms
import math

def MEFromMoments (moms):
    """
    Creates a matrix-exponential distribution that has the
    same moments as given.
    
    Parameters
    ----------
    moms : vector of doubles, length(2*M-1)
        The list of moments. The order of the resulting 
        matrix-exponential distribution is 
        determined based on the number of moments given. To 
        obtain a matrix exponential distribution of order M,
        2*M-1 moments are required.
    
    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial vector of the matrix-exponential 
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential 
        distribution.
    
    References
    ----------
    .. [1] A. van de Liefvoort. The moment problem for 
           continuous distributions. Technical report, 
           University of Missouri, WP-CM-1990-02, Kansas City,
           1990.
    """

    def shift (arr):
        sh = np.roll(arr,1)
        sh[0] = 0
        return sh

    def appie (rmom):
        m = len(rmom)
        if m%2==0:
            rm = rmom[0:m-1]
            m = int (m / 2)
        else:
            rm = rmom
            m = int ((m+1)/2)

        rm = np.insert(np.array(rm),0,1.0)
        f = np.zeros((2*m, 1))
        f[0] = 1.0
        y = np.zeros((2*m, 1))
        n = 0
        k = 0
        q = 1
        d = [0]*m
        alpha = np.zeros((m, m))
        beta = np.zeros((m, 1))
        for i in range(2*m):
            ro = q*np.dot(rm,f)
            nold = n
            n = nold + 1
            yold = y
            if n>0 and ro!=0:
                if k>0:
                    beta[k-1] = ro / np.power(rm[1],d[k-1]+n-1)
                k = k+1
                d[k-1] = n
                n = -n
                q = q / ro
                y = shift(f)
            elif n<=0:
                j = nold + d[k-1]
                alpha[k-1,j] = ro / rm[1]**j
            f = shift(f) - ro*yold

        if sum(d)!=m:
            raise Exception("Insufficient matrix order!")

        K = ml.zeros((m,m))
        K[0,0] = rm[1]
        for i in range(m-1):
            K[i,i+1] = rm[1]
        ind = d[0]
        for i in range (1,m):
            if ind<m:
                inc = d[i]
                ind = ind + inc
                if ind<=m:
                    K[ind-1, ind-inc-d[i-1]] = beta[i-1]
                    for j in range(1,inc+1):
                        K[ind-1, ind-j] = alpha[i,j-1]
        return K

    rmoms = ReducedMomsFromMoms (moms)
    K = appie (rmoms)
    N = int(math.ceil(len(moms)/2.0))
    T = ml.zeros((N,N))
    for i in range(N):
        for j in range(i+1):
            T[i,j] = 1.0
    U = ml.zeros((N,N))
    for i in range(N):
        for j in range(i,N):
            U[i,j] = 1.0 / (N-i)

    alpha = ml.zeros((1,N))
    alpha[0,0] = 1.0
    alpha = alpha * la.inv(T) * U
    A = la.inv(U)*T*K*la.inv(T)*U
    A = la.inv(-A)
    return (alpha, A)

