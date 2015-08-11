butools.dmap.LagCorrelationsFromDRAP
====================================

.. currentmodule:: butools.dmap

.. np:function:: LagCorrelationsFromDRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`acf = LagCorrelationsFromDRAP(H0, H1, L, prec)`
        * - Mathematica:
          - :code:`acf = LagCorrelationsFromDRAP[H0, H1, L, prec]`
        * - Python/Numpy:
          - :code:`acf = LagCorrelationsFromDRAP(H0, H1, L, prec)`

    Returns the lag autocorrelations of a discrete rational 
    arrival process.

    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the discrete rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the discrete rational arrival process
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

    >>> H0 = [0, 0, 0.13; 0, 0.6, 0.18; 0.31, 0.26, 0.02];
    >>> H1 = [0, 1, -0.13; 0, 0.18, 0.04; 0.03, 0.09, 0.29];
    >>> corr = LagCorrelationsFromDRAP(H0,H1,3);
    >>> disp(corr);
         0.014303
        0.0012424
       7.5868e-06

    For Mathematica:

    >>> H0 = {{0, 0, 0.13},{0, 0.6, 0.18},{0.31, 0.26, 0.02}};
    >>> H1 = {{0, 1, -0.13},{0, 0.18, 0.04},{0.03, 0.09, 0.29}};
    >>> corr = LagCorrelationsFromDRAP[H0,H1,3];
    >>> Print[corr];
    {0.01430295723332723, 0.0012424024982963658, 7.586755398724169*^-6}

    For Python/Numpy:

    >>> H0 = ml.matrix([[0, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    >>> H1 = ml.matrix([[0, 1, -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])
    >>> corr = LagCorrelationsFromDRAP(H0,H1,3)
    >>> print(corr)
    [0.014302957233327229, 0.0012424024982963656, 7.5867553987241684e-06]

