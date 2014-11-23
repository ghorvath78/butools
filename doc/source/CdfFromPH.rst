butools.ph.CdfFromPH
====================

.. currentmodule:: butools.ph

.. np:function:: CdfFromPH

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`cdf = CdfFromPH(alpha, A, x, prec)`
        * - Mathematica:
          - :code:`cdf = CdfFromPH[alpha, A, x, prec]`
        * - Python/Numpy:
          - :code:`cdf = CdfFromPH(alpha, A, x, prec)`

    Returns the cummulative distribution function of a
    continuous phase-type distribution.

    Parameters
    ----------
    alpha : matrix, shape (1,M)
        The initial probability vector of the phase-type
        distribution.
    A : matrix, shape (M,M)
        The transient generator matrix of the phase-type
        distribution.
    x : vector of doubles
        The cdf will be computed at these points
    prec : double, optional
        Numerical precision to check if the input PH 
        distribution is valid. The default value is 1e-14.

    Returns
    -------
    cdf : column vector of doubles
        The values of the cdf at the corresponding "x" values

    Examples
    --------    
    For Matlab:
    
    >>> a = [0.1 0.9 0];
    >>> A = [-6.2 2 0; 2 -9 1; 1 0 -3];
    >>> x = (0:0.002:1)';
    >>> cdf = CdfFromPH(a, A, x);
    >>> plot(x, cdf)

    For Python/Numpy:
    
    >>> a=ml.matrix([[0.1, 0.9, 0]])
    >>> A=ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    >>> x = np.linspace(0,1,501)
    >>> cdf = CdfFromPH(a,A,x)
    >>> plt.plot(x,cdf)

