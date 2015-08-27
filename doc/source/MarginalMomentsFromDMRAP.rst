butools.dmap.MarginalMomentsFromDMRAP
=====================================

.. currentmodule:: butools.dmap

.. np:function:: MarginalMomentsFromDMRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`moms = MarginalMomentsFromDMRAP(H, K, precision)`
        * - Mathematica:
          - :code:`moms = MarginalMomentsFromDMRAP[H, K, precision]`
        * - Python/Numpy:
          - :code:`moms = MarginalMomentsFromDMRAP(H, K, precision)`

    Returns the moments of the marginal distribution of a 
    discrete marked rational arrival process.

    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the DMRAP
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

    >>> H0 = [0.15, 0.2, 0.18; -0.23, 0.17, 0.22; 0.19, 0.15, 0.16];
    >>> H1 = [0.01, 0.08, 0.16; 0.02, 0.2, 0.07; 0.02, 0.15, 0.17];
    >>> H2 = [0.14, 0.07, 0.01; 0.19, 0.02, 0.34; 0.06, 0.1, 0];
    >>> moms = MarginalMomentsFromDMRAP({H0, H1, H2});
    >>> disp(moms);
           1.5948       3.4185       9.9595       37.742       177.13

    For Mathematica:

    >>> H0 = {{0.15, 0.2, 0.18},{-0.23, 0.17, 0.22},{0.19, 0.15, 0.16}};
    >>> H1 = {{0.01, 0.08, 0.16},{0.02, 0.2, 0.07},{0.02, 0.15, 0.17}};
    >>> H2 = {{0.14, 0.07, 0.01},{0.19, 0.02, 0.34},{0.06, 0.1, 0}};
    >>> moms = MarginalMomentsFromDMRAP[{H0, H1, H2}];
    >>> Print[moms];
    {1.5947682612697014, 3.418532063993288, 9.959485592524356, 37.74236950164582, 177.12501064870213}

    For Python/Numpy:

    >>> H0 = ml.matrix([[0.15, 0.2, 0.18],[-0.23, 0.17, 0.22],[0.19, 0.15, 0.16]])
    >>> H1 = ml.matrix([[0.01, 0.08, 0.16],[0.02, 0.2, 0.07],[0.02, 0.15, 0.17]])
    >>> H2 = ml.matrix([[0.14, 0.07, 0.01],[0.19, 0.02, 0.34],[0.06, 0.1, 0]])
    >>> moms = MarginalMomentsFromDMRAP([H0, H1, H2])
    >>> print(moms)
    [1.5947682612697014, 3.4185320639932879, 9.9594855925243557, 37.742369501645818, 177.12501064870213]

