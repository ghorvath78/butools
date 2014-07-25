butools.dph.CdfFromMG
=====================

.. currentmodule:: butools.dph

.. np:function:: CdfFromMG

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`cdf = CdfFromMG(alpha, A, x, prec)`
        * - Mathematica:
          - :code:`cdf = CdfFromMG[alpha, A, x, prec]`
        * - Python/Numpy:
          - :code:`cdf = CdfFromMG(alpha, A, x, prec)`

    Returns the cummulative distribution function of a 
    matrix-geometric distribution.

    Parameters
    ----------
    alpha : matrix, shape (1,M)
        The initial vector of the matrix-geometric distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-geometric 
        distribution.
    x : vector of non-negative integers
        The density function will be computed at these points
    prec : double, optional
        Numerical precision to check if the input MG 
        distribution is valid. The default value is 1e-14.

    Returns
    -------
    cdf : column vector of doubles
        The probabilities that the matrix-geometrically 
        distributed random variable is less or equal to
        the corresponding "x" values
        
    Examples
    --------    
    For Matlab:
    
    >>> a=[-0.6 0.3 1.3];
    >>> A=[0.25 0.2 -0.15; 0.3 0.1 0.25; 0 0.2 0.47];
    >>> x = (0:1:20)';
    >>> cdf = CdfFromMG(a, A, x);
    >>> plot(x, cdf)

