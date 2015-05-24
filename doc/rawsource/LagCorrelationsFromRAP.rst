butools.map.LagCorrelationsFromRAP
==================================

.. currentmodule:: butools.map

.. np:function:: LagCorrelationsFromRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`acf = LagCorrelationsFromRAP(H0, H1, L, prec)`
        * - Mathematica:
          - :code:`acf = LagCorrelationsFromRAP[H0, H1, L, prec]`
        * - Python/Numpy:
          - :code:`acf = LagCorrelationsFromRAP(H0, H1, L, prec)`

    Returns the lag autocorrelations of a rational arrival
    process.

    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational arrival process
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
    --------
    For Matlab:
    
    >>> H0=[-2 0 0; 0 -3 1; 0 -1 -2];
    >>> H1=[1.8 0.2 0; 0.2 1.8 0; 0.2 1.8 1];
    >>> LagCorrelationsFromRAP(H0,H1,3)
       -0.0038462
        0.0045604
        0.0058956
    >>> plot(LagCorrelationsFromRAP(H0,H1,20));

    For Python/Numpy:
    
    >>> H0=ml.matrix([[-2, 0, 0],[0, -3, 1],[0, -1, -2]])
    >>> H1=ml.matrix([[1.8, 0.2, 0],[0.2, 1.8, 0],[0.2, 1.8, 1]])
    >>> print(LagCorrelationsFromRAP(H0,H1,3))
    [-0.0038461538461539674, 0.0045604395604395744, 0.005895604395604547]
    >>> plt.plot(LagCorrelationsFromRAP(H0,H1,20))
    
