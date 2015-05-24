butools.dph.PmfFromDPH
======================

.. currentmodule:: butools.dph

.. np:function:: PmfFromDPH

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pmf = PmfFromDPH(alpha, A, x, prec)`
        * - Mathematica:
          - :code:`pmf = PmfFromDPH[alpha, A, x, prec]`
        * - Python/Numpy:
          - :code:`pmf = PmfFromDPH(alpha, A, x, prec)`

    Returns the probability mass function of a discrete
    phase-type distribution.

    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial probability vector of the discrete phase-
        type distribution. The sum of the entries of pi0 is
        less or equal to 1.
    A : matrix, shape (M,M)
        The transient generator matrix of the discrete phase-
        type distribution.
    x : vector of non-negative integers
        The density function will be computed at these points
    prec : double, optional
        Numerical precision to check if the input DPH 
        distribution is valid. The default value is 1e-14.

    Returns
    -------
    pmf : column vector of doubles
        The probabilities that the discrete phase type
        distributed random variable takes the corresponding
        "x" values
        
    Examples
    --------    
    For Matlab:
    
    >>> a=[0.76 0 0.24];
    >>> A=[0.34 0.66 0; 0.79 0.05 0.07; 0.26 0.73 0.01];
    >>> x = (0:1:50)';
    >>> pmf = PmfFromDPH(a, A, x);
    >>> plot(x, pmf)

    For Python/Numpy:
    
    >>> a=ml.matrix([[0.76, 0, 0.24]])
    >>> A=ml.matrix([[0.34, 0.66, 0],[0.79, 0.05, 0.07],[0.26, 0.73, 0.01]])
    >>> x = np.linspace(0,50,51)
    >>> pmf = PmfFromDPH(a,A,x)
    >>> plt.plot(x, pmf)

