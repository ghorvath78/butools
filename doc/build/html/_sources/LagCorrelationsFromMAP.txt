butools.map.LagCorrelationsFromMAP
==================================

.. currentmodule:: butools.map

.. np:function:: LagCorrelationsFromMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`acf = LagCorrelationsFromMAP(D0, D1, L, prec)`
        * - Mathematica:
          - :code:`acf = LagCorrelationsFromMAP[D0, D1, L, prec]`
        * - Python/Numpy:
          - :code:`acf = LagCorrelationsFromMAP(D0, D1, L, prec)`

    Returns the lag autocorrelations of a Markovian arrival
    process.

    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process
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

    >>> D0 = [-5., 0, 1., 1.; 1., -8., 1., 0; 1., 0, -4., 1.; 1., 2., 3., -9.];
    >>> D1 = [0, 1., 0, 2.; 2., 1., 3., 0; 0, 0, 1., 1.; 1., 1., 0, 1.];
    >>> corr = LagCorrelationsFromMAP(D0,D1,3);
    >>> disp(corr);
       0.00012012
       0.00086176
      -0.00022001

    For Mathematica:

    >>> D0 = {{-5., 0, 1., 1.},{1., -8., 1., 0},{1., 0, -4., 1.},{1., 2., 3., -9.}};
    >>> D1 = {{0, 1., 0, 2.},{2., 1., 3., 0},{0, 0, 1., 1.},{1., 1., 0, 1.}};
    >>> corr = LagCorrelationsFromMAP[D0,D1,3];
    >>> Print[corr];
    {0.00012012478025411484, 0.0008617649366101062, -0.00022001393374437001}

    For Python/Numpy:

    >>> D0 = ml.matrix([[-5., 0, 1., 1.],[1., -8., 1., 0],[1., 0, -4., 1.],[1., 2., 3., -9.]])
    >>> D1 = ml.matrix([[0, 1., 0, 2.],[2., 1., 3., 0],[0, 0, 1., 1.],[1., 1., 0, 1.]])
    >>> corr = LagCorrelationsFromMAP(D0,D1,3)
    >>> print(corr)
    [0.00012012478025432307, 0.00086176493661020996, -0.00022001393374405749]

