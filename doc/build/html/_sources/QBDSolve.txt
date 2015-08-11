butools.mam.QBDSolve
====================

.. currentmodule:: butools.mam

.. np:function:: QBDSolve

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[pi0, R] = QBDSolve (B, L, F, L0, prec)`
        * - Mathematica:
          - :code:`{pi0, R} = QBDSolve [B, L, F, L0, prec]`
        * - Python/Numpy:
          - :code:`pi0, R = QBDSolve (B, L, F, L0, prec)`

    Returns the parameters of the matrix-geometrically 
    distributed stationary distribution of a QBD.

    Using vector pi0 and matrix R provided by this function
    the stationary solution can be obtained by

    .. math::
        \pi_k=\pi_0 R^k.    
    
    Parameters
    ----------
    B : matrix, shape (N,N)
        The matrix corresponding to backward transitions
    L : matrix, shape (N,N)
        The matrix corresponding to local transitions
    F : matrix, shape (N,N)
        The matrix corresponding to forward transitions
    L0 : matrix, shape (N,N)
        The matrix corresponding to local transitions at
        level zero
    precision : double, optional
        The fundamental matrix R is computed up to this
        precision. The default value is 1e-14
    
    Returns
    -------
    pi0 : matrix, shape (1,N)
        The stationary probability vector of level zero
    R : matrix, shape (N,N)
        The matrix parameter of the matrix geometrical
        distribution of the QBD 
    
    Examples
    --------
    For Matlab:
    
    >>> B = [0,0;3,4];
    >>> L = [-6,5;3,-12];
    >>> F = [1,0;2,0];
    >>> L0 = [-6,5;6,-8];
    >>> [pi0, R] = QBDSolve (B, L, F, L0);
    >>> pi0
          0.22992      0.18681
    >>> R
          0.27839      0.14286
          0.55678      0.28571
    >>> pi1 = pi0*R
          0.16802     0.086221
    >>> pi2 = pi0*R^2
         0.094781     0.048638

    For Python/Numpy:
    
    >>> B = ml.matrix([[0,0],[3,4]])
    >>> L = ml.matrix([[-6,5],[3,-12]])
    >>> F = ml.matrix([[1,0],[2,0]])
    >>> L0 = ml.matrix([[-6,5],[6,-8]])
    >>> pi0, R = QBDSolve (B, L, F, L0)
    >>> print(pi0)
    [[ 0.22992392  0.18681319]]
    >>> print(R)
    [[ 0.27838828  0.14285714]
     [ 0.55677656  0.28571429]]
    >>> print(pi0*R)
    [[ 0.16802133  0.08622147]]
    >>> print(pi0*R**2)
    [[ 0.09478126  0.04863775]]

