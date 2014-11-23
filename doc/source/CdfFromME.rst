butools.ph.CdfFromME
====================

.. currentmodule:: butools.ph

.. np:function:: CdfFromME

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`cdf = CdfFromME(alpha, A, x, prec)`
        * - Mathematica:
          - :code:`cdf = CdfFromME[alpha, A, x, prec]`
        * - Python/Numpy:
          - :code:`cdf = CdfFromME(alpha, A, x, prec)`

    Returns the cummulative distribution function of a
    matrix-exponential distribution.

    Parameters
    ----------
    alpha : matrix, shape (1,M)
        The initial vector of the matrix-exponential
        distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential
        distribution.
    x : vector of doubles
        The cdf will be computed at these points
    prec : double, optional
        Numerical precision to check if the input ME 
        distribution is valid. The default value is 1e-14.

    Returns
    -------
    cdf : column vector of doubles
        The values of the cdf at the corresponding "x" values

    Examples
    --------    
    For Matlab:
    
    >>> a = [0.2, 0.3, 0.5];
    >>> A = [-1,0,0;0,-3,2;0,-2,-3];
    >>> x = (0:0.01:2)';
    >>> cdf = CdfFromME(a, A, x);
    >>> plot(x, cdf)

    For Python/Numpy:
    
    >>> a = ml.matrix([[0.2, 0.3, 0.5]])
    >>> A = ml.matrix([[-1,0,0],[0,-3,2],[0,-2,-3]])
    >>> x = np.linspace(0,2,201)
    >>> cdf = CdfFromME(a,A,x)
    >>> plt.plot(x,cdf)

