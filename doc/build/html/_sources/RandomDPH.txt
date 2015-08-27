butools.dph.RandomDPH
=====================

.. currentmodule:: butools.dph

.. np:function:: RandomDPH

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = RandomDPH(order, mean, zeroEntries, maxTrials, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = RandomDPH[order, mean, zeroEntries, maxTrials, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = RandomDPH(order, mean, zeroEntries, maxTrials, prec)`

    Returns a random discrete phase-type distribution with a 
    given mean value.

    Parameters
    ----------
    order : int
        The size of the discrete phase-type distribution
    mean : double, optional
        The mean of the discrete phase-type distribution 
    zeroEntries : int, optional
        The number of zero entries in the initial vector, 
        generator matrix and closing vector
    maxTrials : int, optional
        The maximum number of trials to find a proper DPH 
        (that has an irreducible phase process and none of 
        its parameters is all-zero). The default value is 
        1000.
    prec : double, optional
        Numerical precision for checking the irreducibility.
        The default value is 1e-14.

    Returns
    -------
    alpha : vector, shape (1,M)
        The initial probability vector of the phase-type 
        distribution.
    A : matrix, shape (M,M)
        The transient generator matrix of the phase-type 
        distribution.

    Notes
    -----
    If the procedure fails, try to increase the 'maxTrials'
    parameter, or increase the mean value.

    Examples
    ========
    For Matlab:

    >>> [a, A] = RandomDPH(3, 10, 5);
    >>> disp(a);
          0.39611      0.38982      0.21407
    >>> disp(A);
           0.8489     0.031291            0
          0.59994      0.40006            0
          0.37646            0      0.62354

    For Mathematica:

    >>> {a, A} = RandomDPH[3, 10, 5];
    "CheckProbMatrix: the matrix has negative element (at precision "1.*^-7")!"
    >>> Print[a];
    {0, 0, 1}
    >>> Print[A];
    {{0.8271227985086882, 0.06664258193377719, 0.},
     {0., 0.9130962745792692, 0.},
     {0.05885243293901752, 0.02706649617956366, 0.8087057205828357}}

    For Python/Numpy:

    >>> a, A = RandomDPH(3, 10, 5)
    >>> print(a)
    [[ 1.  0.  0.]]
    >>> print(A)
    [[ 0.67742  0.09221  0.07805]
     [ 0.08848  0.83115  0.     ]
     [ 0.1876   0.       0.8124 ]]

