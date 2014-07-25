butools.map.MarginalMomentsFromMMAP
===================================

.. currentmodule:: butools.map

.. np:function:: MarginalMomentsFromMMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`moms = MarginalMomentsFromMMAP(D, K, precision)`
        * - Mathematica:
          - :code:`moms = MarginalMomentsFromMMAP[D, K, precision]`
        * - Python/Numpy:
          - :code:`moms = MarginalMomentsFromMMAP(D, K, precision)`

    Returns the moments of the marginal distribution of a 
    marked Markovian arrival process.

    Parameters
    ----------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the MMAP
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
    
    >>> D0=[-1.78 0.29; 0.07 -0.92];
    >>> D1=[0.15 0.49; 0.23 0.36];
    >>> D2=[0.11 0.2; 0.01 0];
    >>> D3=[0.14 0.4; 0.11 0.14];
    >>> MarginalMomentsFromMMAP({D0,D1,D2,D3})
           1.0007       2.1045       6.8277

