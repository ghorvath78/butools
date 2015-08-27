butools.dmap.LagCorrelationsFromDMAP
====================================

.. currentmodule:: butools.dmap

.. np:function:: LagCorrelationsFromDMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`acf = LagCorrelationsFromDMAP(D0, D1, L, prec)`
        * - Mathematica:
          - :code:`acf = LagCorrelationsFromDMAP[D0, D1, L, prec]`
        * - Python/Numpy:
          - :code:`acf = LagCorrelationsFromDMAP(D0, D1, L, prec)`

    Returns the lag autocorrelations of a discrete Markovian 
    arrival process.

    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the discrete Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the discrete Markovian arrival process
    L : double, optional
        The number of lags to compute. The default value is 1
    prec : double, optional
        Numerical precision to check if the input is valid. 
        The default value is 1e-14

    Returns
    -------
    acf : column vector of doubles, length (L)
        The lag autocorrelation function up to lag L
        
    Examples
    ========
    For Matlab:

    >>> D0 = [0, 0.02, 0, 0; 0, 0.17, 0.2, 0.14; 0.16, 0.17, 0.02, 0.18; 0, 0, 0, 0.12];
    >>> D1 = [0, 0.88, 0.1, 0; 0.18, 0.07, 0.14, 0.1; 0.13, 0.15, 0.15, 0.04; 0.31, 0.18, 0.12, 0.27];
    >>> corr = LagCorrelationsFromDMAP(D0, D1, 3);
    >>> disp(corr);
        -0.045859     0.010753   -0.0047996

    For Mathematica:

    >>> D0 = {{0, 0.02, 0, 0},{0, 0.17, 0.2, 0.14},{0.16, 0.17, 0.02, 0.18},{0, 0, 0, 0.12}};
    >>> D1 = {{0, 0.88, 0.1, 0},{0.18, 0.07, 0.14, 0.1},{0.13, 0.15, 0.15, 0.04},{0.31, 0.18, 0.12, 0.27}};
    >>> corr = LagCorrelationsFromDMAP[D0, D1, 3];
    >>> Print[corr];
    {-0.04585887310401268, 0.010753286512163932, -0.00479959597519405}

    For Python/Numpy:

    >>> D0 = ml.matrix([[0, 0.02, 0, 0],[0, 0.17, 0.2, 0.14],[0.16, 0.17, 0.02, 0.18],[0, 0, 0, 0.12]])
    >>> D1 = ml.matrix([[0, 0.88, 0.1, 0],[0.18, 0.07, 0.14, 0.1],[0.13, 0.15, 0.15, 0.04],[0.31, 0.18, 0.12, 0.27]])
    >>> corr = LagCorrelationsFromDMAP(D0, D1, 3)
    >>> print(corr)
    [-0.04586  0.01075 -0.0048 ]

