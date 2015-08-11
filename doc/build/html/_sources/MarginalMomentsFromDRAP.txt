butools.dmap.MarginalMomentsFromDRAP
====================================

.. currentmodule:: butools.dmap

.. np:function:: MarginalMomentsFromDRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`moms = MarginalMomentsFromDRAP(H0, H1, K, precision)`
        * - Mathematica:
          - :code:`moms = MarginalMomentsFromDRAP[H0, H1, K, precision]`
        * - Python/Numpy:
          - :code:`moms = MarginalMomentsFromDRAP(H0, H1, K, precision)`

    Returns the moments of the marginal distribution of a 
    discrete rational arrival process.

    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the discrete rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the discrete rational arrival process
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

    >>> H0 = [0, 0, 0.13; 0, 0.6, 0.18; 0.31, 0.26, 0.02];
    >>> H1 = [0, 1, -0.13; 0, 0.18, 0.04; 0.03, 0.09, 0.29];
    >>> moms = MarginalMomentsFromDRAP(H0,H1);
    >>> disp(moms);
            3.207       16.898       130.77       1347.1        17343

    For Mathematica:

    >>> H0 = {{0, 0, 0.13},{0, 0.6, 0.18},{0.31, 0.26, 0.02}};
    >>> H1 = {{0, 1, -0.13},{0, 0.18, 0.04},{0.03, 0.09, 0.29}};
    >>> moms = MarginalMomentsFromDRAP[H0,H1];
    >>> Print[moms];
    {3.20702366840782, 16.897636691953394, 130.7705457435602, 1347.0743893905096, 17343.182467560622}

    For Python/Numpy:

    >>> H0 = ml.matrix([[0, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    >>> H1 = ml.matrix([[0, 1, -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])
    >>> moms = MarginalMomentsFromDRAP(H0,H1)
    >>> print(moms)
    [3.2070236684078202, 16.897636691953394, 130.77054574356021, 1347.0743893905096, 17343.182467560622]

