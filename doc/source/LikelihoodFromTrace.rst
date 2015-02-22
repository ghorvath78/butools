butools.fitting.LikelihoodFromTrace
===================================

.. currentmodule:: butools.fitting

.. np:function:: LikelihoodFromTrace

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`logli = LikelihoodFromTrace(trace, X, Y, prec)`
        * - Mathematica:
          - :code:`logli = LikelihoodFromTrace[trace, X, Y, prec]`
        * - Python/Numpy:
          - :code:`logli = LikelihoodFromTrace(trace, X, Y, prec)`

    Evaluates the log-likelihood of a trace with the given PH
    distribution or MAP. The result is divided by the length
    of the trace.
    
    If X is a row vector, than (X,Y) is interpreted as a PH
    distribution, otherwise (X,Y) is considered to be a MAP.
    
    Parameters
    ----------
    trace : column vector, length K
        The samples of the trace
    X : matrix, shape (1,M) or (M,M)
        If X is a row vector, it is the initial probability
        vector of the PH distribution. If X is a square
        matrix, it is interpreted as the D0 matrix of a MAP
    Y : matrix, (M,M)
        If X is a row vector, Y is the transient generator
        of the PH distribution. If X is a square matrix, Y
        is interpreted as the D1 matrix of a MAP
    prec : double, optional
        Numerical precision used by the randomization. The
        default value is 1e-14.

    Returns
    -------
    logli : double
        The log likelihood divided by the size of the trace
        
    Notes
    -----
    The procedure is much faster with PH distributions.

    Examples
    --------    
    For Matlab:
    
    >>> tr = dlmread('trace.txt');
    >>> moms = MarginalMomentsFromTrace(tr,3);
    >>> [alpha,A] = APHFrom3Moments(moms);
    >>> LikelihoodFromTrace(tr,alpha,A)
           4.8496
    >>> corr1 = LagCorrelationsFromTrace(tr,1);
    >>> [D0,D1]=MAPFromFewMomentsAndCorrelations(moms,corr1)
    >>> LikelihoodFromTrace(tr,D0,D1)
           4.6523    
    
    For Python/Numpy:
    
    >>> tr = np.loadtxt('trace.txt')
    >>> moms = MarginalMomentsFromTrace(tr,3)
    >>> alpha,A = APHFrom3Moments(moms)
    >>> print(LikelihoodFromTrace(tr,alpha,A))
    4.8496353019637972
    >>> corr1 = LagCorrelationsFromTrace(tr,1)[0]
    >>> D0,D1=MAPFromFewMomentsAndCorrelations(moms,corr1)
    >>> LikelihoodFromTrace(tr,D0,D1)
    4.652346100436191

           
    
