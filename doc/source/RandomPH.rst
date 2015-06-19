butools.ph.RandomPH
===================

.. currentmodule:: butools.ph

.. np:function:: RandomPH

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = RandomPH(order, mean, zeroEntries, maxTrials, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = RandomPH[order, mean, zeroEntries, maxTrials, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = RandomPH(order, mean, zeroEntries, maxTrials, prec)`

    Returns a random phase-type distribution with a given 
    order.

    Parameters
    ----------
    order : int
        The size of the phase-type distribution
    mean : double, optional
        The mean of the phase-type distribution 
    zeroEntries : int, optional
        The number of zero entries in the initial vector, 
        generator matrix and closing vector
    maxTrials : int, optional
        The maximum number of trials to find a proper PH 
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
    parameter.   

    Examples
    ========
    For Matlab:

    >>> [a,A] = RandomPH(3,8,4);
    >>> disp(a);
          0.28793      0.34404      0.36802
    >>> disp(A);
           -2.115      0.40925       1.1338
           1.9391      -1.9391            0
                0      0.52322     -0.52322

    For Mathematica:

    >>> {a,A} = RandomPH[3,8,4];
    >>> Print[a];
    {0, 1, 0}
    >>> Print[A];
    {{-0.4376516294910215, 0.05516364914230338, 0.16082329157888142},
     {0.17964205281519655, -0.2942958850413159, 0.},
     {0.13943848467964445, 0.12380908822485874, -0.2632475729045032}}

    For Python/Numpy:

    >>> a,A = RandomPH(3,8,4)
    >>> print(a)
    [[ 0.53323  0.36287  0.10391]]
    >>> print(A)
    [[-1.79031  0.85925  0.34048]
     [ 0.1859  -0.1859   0.     ]
     [ 0.       0.      -0.67749]]

