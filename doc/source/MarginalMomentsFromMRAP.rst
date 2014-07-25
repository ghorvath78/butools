butools.map.MarginalMomentsFromMRAP
===================================

.. currentmodule:: butools.map

.. np:function:: MarginalMomentsFromMRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`moms = MarginalMomentsFromMRAP(H, K, precision)`
        * - Mathematica:
          - :code:`moms = MarginalMomentsFromMRAP[H, K, precision]`
        * - Python/Numpy:
          - :code:`moms = MarginalMomentsFromMRAP(H, K, precision)`

    Returns the moments of the marginal distribution of a 
    marked rational arrival process.

    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the MRAP
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
    
    >>> H0=[-5 0.28 0.9 1; 1 -8 0.9 0.1; 0.9 0.1 -4 1; 1 2 3 -9]
    >>> H1=[-0.08 0.7 0.1 0.1; 0.1 1 1.8 0.1; 0.1 0.1 0.1 0.7; 0.7 0.1 0.1 0.1]
    >>> H2=[0.1 0.1 0.1 1.7; 1.8 0.1 1 0.1; 0.1 0.1 0.7 0.1; 0.1 1 0.1 0.8]
    >>> MarginalMomentsFromMRAP({H0,H1,H2})
          0.33951      0.24583      0.27424      0.41206      0.77677       1.7594       4.6515

