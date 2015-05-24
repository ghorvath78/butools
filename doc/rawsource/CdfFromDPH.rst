butools.dph.CdfFromDPH
======================

.. currentmodule:: butools.dph

.. np:function:: CdfFromDPH

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`cdf = CdfFromDPH(alpha, A, x, prec)`
        * - Mathematica:
          - :code:`cdf = CdfFromDPH[alpha, A, x, prec]`
        * - Python/Numpy:
          - :code:`cdf = CdfFromDPH(alpha, A, x, prec)`

    Returns the cummulative distribution function of a 
    discrete phase-type distribution.

    Parameters
    ----------
    alpha : matrix, shape (1,M)
        The initial probability vector of the discrete phase-
        type distribution.
    A : matrix, shape (M,M)
        The transition probability  matrix of the discrete phase-
        type distribution.
    x : vector of non-negative integers
        The cdf will be computed at these points
    prec : double, optional
        Numerical precision to check if the input DPH 
        distribution is valid. The default value is 1e-14.

    Returns
    -------
    cdf : column vector of doubles
        The probabilities that the discrete phase type 
        distributed random variable is less or equal to the
        corresponding "x" values
        
    Examples
    --------    
    For Matlab:
    
    >>> a=[0.76 0 0.24];
    >>> A=[0.34 0.66 0; 0.79 0.05 0.07; 0.26 0.73 0.01];
    >>> x = (0:1:100)';
    >>> cdf = CdfFromDPH(a, A, x);
    >>> plot(x, cdf)

    For Python/Numpy:
    
    >>> a=ml.matrix([[0.76, 0, 0.24]])
    >>> A=ml.matrix([[0.34, 0.66, 0],[0.79, 0.05, 0.07],[0.26, 0.73, 0.01]])
    >>> x = np.linspace(0,100,101)
    >>> pmf = CdfFromDPH(a,A,x)
    >>> plt.plot(x, pmf)

