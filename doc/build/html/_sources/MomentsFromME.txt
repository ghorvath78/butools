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
    ========
    For Matlab:

    >>> a = [0.2, 0.3, 0.5];
    >>> A = [-1, 0, 0; 0, -3, 2; 0, -2, -3];
    >>> moms = MomentsFromME(a,A);
    >>> disp(moms);
          0.35385      0.41893       1.1552       4.6998       23.838
    >>> moms = MomentsFromME(a,A,9);
    >>> disp(moms);
          0.35385      0.41893       1.1552       4.6998       23.838       143.78       1007.8       8064.3        72578

    For Mathematica:

    >>> a = {0.2, 0.3, 0.5};
    >>> A = {{-1, 0, 0},{0, -3, 2},{0, -2, -3}};
    >>> moms = MomentsFromME[a,A];
    >>> Print[moms];
    {0.35384615384615387, 0.41893491124260357, 1.1552116522530724, 4.699835439935577, 23.837756165615836}
    >>> moms = MomentsFromME[a,A,9];
    >>> Print[moms];
    {0.35384615384615387, 0.41893491124260357, 1.1552116522530724, 4.699835439935577, 23.837756165615836, 143.78185836646944, 1007.8194071104502, 8064.272882521486, 72578.13371878522}

    For Python/Numpy:

    >>> a = ml.matrix([[0.2, 0.3, 0.5]])
    >>> A = ml.matrix([[-1, 0, 0],[0, -3, 2],[0, -2, -3]])
    >>> moms = MomentsFromME(a,A)
    >>> print(moms)
    [0.35384615384615381, 0.41893491124260357, 1.1552116522530724, 4.6998354399355771, 23.837756165615836]
    >>> moms = MomentsFromME(a,A,9)
    >>> print(moms)
    [0.35384615384615381, 0.41893491124260357, 1.1552116522530724, 4.6998354399355771, 23.837756165615836, 143.78185836646944, 1007.8194071104502, 8064.2728825214863, 72578.133718785219]

