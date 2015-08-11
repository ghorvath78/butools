butools.map.MarginalMomentsFromRAP
==================================

.. currentmodule:: butools.map

.. np:function:: MarginalMomentsFromRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`moms = MarginalMomentsFromRAP(H0, H1, K, precision)`
        * - Mathematica:
          - :code:`moms = MarginalMomentsFromRAP[H0, H1, K, precision]`
        * - Python/Numpy:
          - :code:`moms = MarginalMomentsFromRAP(H0, H1, K, precision)`

    Returns the moments of the marginal distribution of a 
    rational arrival process.

    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational arrival process
    K : int, optional
        Number of moments to compute. If K=0, 2*M-1 moments
        are computed. The default value is K=0.
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14

    Returns
    -------
    moms : row vector of doubles, length K
        The vector of moments.

    Examples
    ========
    For Matlab:

    >>> H0 = [-2., 0, 0; 0, -3., 1.; 0, -1., -2.];
    >>> H1 = [1.8, 0.2, 0; 0.2, 1.8, 0; 0.2, 1.8, 1.];
    >>> moms = MarginalMomentsFromRAP(H0,H1);
    >>> disp(moms);
          0.44444      0.38095      0.48299      0.82216       1.7944

    For Mathematica:

    >>> H0 = {{-2., 0, 0},{0, -3., 1.},{0, -1., -2.}};
    >>> H1 = {{1.8, 0.2, 0},{0.2, 1.8, 0},{0.2, 1.8, 1.}};
    >>> moms = MarginalMomentsFromRAP[H0,H1];
    >>> Print[moms];
    {0.4444444444444444, 0.380952380952381, 0.48299319727891166, 0.8221574344023325, 1.794391225878107}

    For Python/Numpy:

    >>> H0 = ml.matrix([[-2., 0, 0],[0, -3., 1.],[0, -1., -2.]])
    >>> H1 = ml.matrix([[1.8, 0.2, 0],[0.2, 1.8, 0],[0.2, 1.8, 1.]])
    >>> moms = MarginalMomentsFromRAP(H0,H1)
    >>> print(moms)
    [0.44444444444444442, 0.38095238095238093, 0.48299319727891149, 0.82215743440233213, 1.7943912258781058]

