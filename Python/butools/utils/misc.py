import numpy as np
import numpy.matlib as ml
import numpy.linalg as la

def Vec(A):
    """
    Column stacking (vec operator).

    Parameters
    ----------
    A : matrix, shape (M,M)

    Returns
    -------
    v : matrix, shape (M,1)
        v constists of the columns of A stacked under each other
    """
    return np.reshape(A.T, (-1,1))
   
def Diag(v):
    """
    This function works with vectors and matrices as well.
    
    In case of square matrices:

    Parameters
    ----------
    v : matrix, shape (M,M) or (1,M) or (M,1)
    
    Returns
    -------
    d : matrix, shape (M,1) of (M,M)
        If v is a square matrix, d is a column vector of the diagonal elements of matrix v.
        If v is a row or a column vector, d is a diagonal matrix constructed from the elements of v.
    """
    if len(v.shape)>1 and v.shape[0]==v.shape[1]:
        return ml.matrix(np.diag(v)).T
    elif len(v.shape)==1:
        return ml.matrix(np.diag(v))
    else:
        return ml.matrix(np.diag(v.A.flatten()))

def Linsolve(A,b):
    """
    Solves the linear system A*x=b (if b is a column vector), or x*A=b (if b is 
    a row vector).
    
    Matrix "A" does not need to be square, this function uses rank-revealing
    QR decomposition to solve the system.
    
    Parameters
    ----------
    A : matrix, shape (M,N)
        The coefficient matrix of the linear system.
    b : matrix, shape (M,1) or (1,N)
        The right hand side of the linear system
        
    Returns
    -------
    x : matrix, shape (M,1) or (1,N)
        If b is a column vector, then x is the solution of A*x=b.       
        If b is a row vector, it returns the solution of x*A=b.
    """
    if b.shape[0]==1:
        x = Linsolve(np.conj(A.T), np.conj(b.T))
        return np.conj(x.T)
    elif b.shape[1]==1:
        Q,R = la.qr(A)
        N = A.shape[1]
        return ml.matrix(la.solve(R[0:N,0:N], np.array(np.conj(Q.T)*b).flatten()[0:N])).T

def SumMatrixList(C):
    sumC = ml.zeros(C[0].shape)
    for i in range(len(C)):
        sumC += C[i]
    return sumC

def Length(v):
    if isinstance(v,np.ndarray):
        if v.shape[0]==1:
            return v.shape[1]
        else:
            return v.shape[0]
    else:
        return len(v)
