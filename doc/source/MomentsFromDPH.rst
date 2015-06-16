butools.dph.MomentsFromDPH
==========================

.. currentmodule:: butools.dph

.. np:function:: MomentsFromDPH

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`moms = MomentsFromDPH(alpha, A, K, prec)`
        * - Mathematica:
          - :code:`moms = MomentsFromDPH[alpha, A, K, prec]`
        * - Python/Numpy:
          - :code:`moms = MomentsFromDPH(alpha, A, K, prec)`

    Returns the first K moments of a discrete phase-type
    distribution.

    Parameters
    ----------
    alpha : vector, shape (1,M)
        The initial probability vector of the discrete phase-
        type distribution. The sum of the entries of pi0 is 
        less or equal to 1.
    A : matrix, shape (M,M)
        The transient generator matrix of the discrete phase-
        type distribution.
    K : int, optional
        Number of moments to compute. If K=0, 2*M-1 moments
        are computed. The default value is 0.
    prec : double, optional
        Numerical precision for checking the input.
        The default value is 1e-14.

    Returns
    -------
    moms : row vector of doubles
        The vector of moments.
        
    Examples
    ========
    For Matlab:

    >>> a = [0.76, 0, 0.24];
    >>> A = [0.34, 0.66, 0; 0.79, 0.05, 0.07; 0.26, 0.73, 0.01];
    >>> moms = MomentsFromDPH(a,A,5);
    >>> disp(moms);
           26.995         1398   1.0853e+05   1.1233e+07   1.4533e+09

    For Mathematica:

    >>> a = {0.76, 0, 0.24};
    >>> A = {{0.34, 0.66, 0},{0.79, 0.05, 0.07},{0.26, 0.73, 0.01}};
    >>> moms = MomentsFromDPH[a,A,5];
    >>> Print[moms];
    {26.995340611502304, 1397.9993695881547, 108525.47866809377, 1.1232963460675944*^7, 1.4533393399621515*^9}

    For Python/Numpy:

    >>> a = ml.matrix([[0.76, 0, 0.24]])
    >>> A = ml.matrix([[0.34, 0.66, 0],[0.79, 0.05, 0.07],[0.26, 0.73, 0.01]])
    >>> moms = MomentsFromDPH(a,A,5)
    >>> print(moms)
    [26.995340611502307, 1397.9993695881547, 108525.47866809377, 11232963.460675946, 1453339339.9621518]

