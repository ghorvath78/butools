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
    ========
    For Matlab:

    >>> B = [0., 0.; 3., 4.];
    >>> L = [-6., 5.; 3., -12.];
    >>> F = [1., 0.; 2., 0.];
    >>> L0 = [-6., 5.; 6., -8.];
    >>> [pi0, R] = QBDSolve(B, L, F, L0);
    >>> disp(pi0);
          0.22992      0.18681
    >>> disp(R);
          0.27839      0.14286
          0.55678      0.28571

    For Mathematica:

    
    For Python/Numpy:

    >>> B = ml.matrix([[0., 0.],[3., 4.]])
    >>> L = ml.matrix([[-6., 5.],[3., -12.]])
    >>> F = ml.matrix([[1., 0.],[2., 0.]])
    >>> L0 = ml.matrix([[-6., 5.],[6., -8.]])
    >>> pi0, R = QBDSolve(B, L, F, L0)
    Final Residual Error for G:  1.38777878078e-16
    Final Residual Error for R:  5.55111512313e-17
    >>> print(pi0)
    [[ 0.22992  0.18681]]
    >>> print(R)
    [[ 0.27839  0.14286]
     [ 0.55678  0.28571]]

