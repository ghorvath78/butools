butools.mam.QBDSolve
====================

.. currentmodule:: butools.mam

.. np:function:: QBDSolve

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[pi0, R] = QBDSolve (B, L, F, L0)`
        * - Mathematica:
          - :code:`{pi0, R} = QBDSolve [B, L, F, L0]`
        * - Python/Numpy:
          - :code:`pi0, R = QBDSolve (B, L, F, L0)`

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
    L : matrix, shape (N,N)
        The matrix corresponding to local transitions at
        level zero
    
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

