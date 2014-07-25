butools.map.CanonicalFromMAP2
=============================

.. currentmodule:: butools.map

.. np:function:: CanonicalFromMAP2

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[G0, G1] = CanonicalFromMAP2(D0, D1, prec)`
        * - Mathematica:
          - :code:`{G0, G1} = CanonicalFromMAP2[D0, D1, prec]`
        * - Python/Numpy:
          - :code:`G0, G1 = CanonicalFromMAP2(D0, D1, prec)`

    Returns the canonical form of an order-2 Markovian
    arrival process.

    Parameters
    ----------
    D0 : matrix, shape (2,2)
        The D0 matrix of the MAP(2)
    D1 : matrix, shape (2,2)
        The D1 matrix of the MAP(2)
    prec : double, optional
        Numerical precision to check the input, default 
        value is 1e-14

    Returns
    -------
    G0 : matrix, shape (1,2)
        The D0 matrix of the canonical MAP(2)
    G1 : matrix, shape (2,2)
        The D1 matrix of the canonical MAP(2)

    Notes
    -----
    This procedure calculates 3 marginal moments and the lag-1
    autocorrelation of the input and calls 'MAP2FromMoments'.

    Examples
    --------
    For Matlab:
    
    >>> D0=[-14 1; 1 -25];
    >>> D1=[6 7; 3 21];
    >>> [H0,H1]=CanonicalFromMAP2(D0,D1);
    >>> H0
           -13.91        9.199
                0       -25.09
    >>> H1
           4.7108            0
            2.801       22.289
    >>> C=SimilarityMatrix(H0,D0);
    >>> dissimilarity = norm(H0*C-C*D0) + norm(H1*C-C*D1)
          5.3e-13
  

