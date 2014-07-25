butools.ph.MomentsFromME
========================

.. currentmodule:: butools.ph

.. np:function:: MomentsFromME

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`moms = MomentsFromME(alpha, A, K, prec)`
        * - Mathematica:
          - :code:`moms = MomentsFromME[alpha, A, K, prec]`
        * - Python/Numpy:
          - :code:`moms = MomentsFromME(alpha, A, K, prec)`

    Returns the first K moments of a matrix-exponential
    distribution.

    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial vector of the matrix-exponential
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential
        distribution.
    K : int, optional
        Number of moments to compute. If K=0, 2*M-1 moments
        are computed. The default value is K=0.
    prec : double, optional
        Numerical precision for checking the input.
        The default value is 1e-14.
    
    Returns
    -------
    moms : row vector of doubles
        The vector of moments.
        
    Examples
    --------    
    For Matlab:
    
    >>> a = [0.2, 0.3, 0.5];
    >>> A = [-1,0,0;0,-3,2;0,-2,-3];
    >>> moms = MomentsFromME(a,A)
          0.35385      0.41893       1.1552       4.6998       23.838
    >>> mean = MomentsFromME(a,A,1)
          0.35385

