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
    ========
    For Matlab:

    >>> tr = dlmread('/home/gabor/github/butools/test/data/bctrace.iat');
    >>> [alpha, A] = APHFrom3Moments(MarginalMomentsFromTrace(tr, 3));
    >>> [D0, D1] = MAPFromFewMomentsAndCorrelations(MarginalMomentsFromTrace(tr, 3), LagCorrelationsFromTrace(tr, 1));
    >>> logliPH = LikelihoodFromTrace(tr, alpha, A);
    >>> disp(logliPH);
           4.8496
    >>> logliMAP = LikelihoodFromTrace(tr, D0, D1);
    >>> disp(logliMAP);
           4.6523

    For Mathematica:

    >>> tr = Flatten[Import["/home/gabor/github/butools/test/data/bctrace.iat","CSV"]];
    >>> {alpha, A} = APHFrom3Moments[MarginalMomentsFromTrace[tr, 3]];
    >>> {D0, D1} = MAPFromFewMomentsAndCorrelations[MarginalMomentsFromTrace[tr, 3], LagCorrelationsFromTrace[tr, 1][[1]]];
    >>> logliPH = LikelihoodFromTrace[tr, alpha, A];
    >>> Print[logliPH];
    4.849635301970579
    >>> logliMAP = LikelihoodFromTrace[tr, D0, D1];
    >>> Print[logliMAP];
    4.65234643278208

    For Python/Numpy:

    >>> tr = np.loadtxt("/home/gabor/github/butools/test/data/bctrace.iat")
    >>> alpha, A = APHFrom3Moments(MarginalMomentsFromTrace(tr, 3))
    >>> D0, D1 = MAPFromFewMomentsAndCorrelations(MarginalMomentsFromTrace(tr, 3), LagCorrelationsFromTrace(tr, 1)[0])
    >>> logliPH = LikelihoodFromTrace(tr, alpha, A)
    >>> print(logliPH)
    4.84963530196
    >>> logliMAP = LikelihoodFromTrace(tr, D0, D1)
    >>> print(logliMAP)
    4.652346100436191

