butools.ph.MomentsFromPH
========================

.. currentmodule:: butools.ph

.. np:function:: MomentsFromPH

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`moms = MomentsFromPH(alpha, A, K, prec)`
        * - Mathematica:
          - :code:`moms = MomentsFromPH[alpha, A, K, prec]`
        * - Python/Numpy:
          - :code:`moms = MomentsFromPH(alpha, A, K, prec)`

    Returns the first K moments of a continuous phase-type
    distribution.

    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial probability vector of the phase-type
        distribution.
    A : matrix, shape (M,M)
        The transient generator matrix of the phase-type
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
    ========
    For Matlab:

    >>> a = [0.1,0.9,0];
    >>> A = [-6.2, 2, 0; 2, -9, 1; 1, 0, -3];
    >>> moms = MomentsFromPH(a, A, 5);
    >>> disp(moms);
          0.20939      0.10449     0.089092      0.11027      0.17953

    For Mathematica:

    
    For Python/Numpy:

    >>> a = ml.matrix([[0.1,0.9,0]])
    >>> A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    >>> moms = MomentsFromPH(a, A, 5)
    >>> print(moms)
    [0.20938722294654497, 0.10448912014333091, 0.089091500391907288, 0.11026674096545433, 0.17953027324720897]

