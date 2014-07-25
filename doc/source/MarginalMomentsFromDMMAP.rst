butools.dmap.MarginalMomentsFromDMMAP
=====================================

.. currentmodule:: butools.dmap

.. np:function:: MarginalMomentsFromDMMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`moms = MarginalMomentsFromDMMAP(D, K, precision)`
        * - Mathematica:
          - :code:`moms = MarginalMomentsFromDMMAP[D, K, precision]`
        * - Python/Numpy:
          - :code:`moms = MarginalMomentsFromDMMAP(D, K, precision)`

    Returns the moments of the marginal distribution of a 
    discrete marked Markovian arrival process.

    Parameters
    ----------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the DMMAP
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
    --------
    For Matlab:

    >>> D0=[0.34 0 0; 0.06 0.05 0.03; 0.11 0.13 0];
    >>> D1=[0.3 0 0; 0.16 0.18 0.05; 0.15 0.04 0.09];
    >>> D2=[0 0.01 0; 0.1 0.07 0.08; 0.13 0.12 0.13];
    >>> D3=[0.35 0 0; 0 0.18 0.04; 0.06 0.03 0.01];
    >>> MarginalMomentsFromDMMAP({D0,D1,D2,D3})
           1.5037       3.0278       8.4243       31.097       143.88

