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

    >>> [a, A] = RandomPH(3, 8, 4);
    >>> disp(a);
          0.50814            0      0.49186
    >>> disp(A);
         -0.38196      0.17775      0.17535
                0     -0.24494       0.1677
                0            0      -0.1555

    For Mathematica:

    >>> {a, A} = RandomPH[3, 8, 4];
    >>> Print[a];
    {0.31004439863858146, 0.047613750639457386, 0.6423418507219611}
    >>> Print[A];
    {{-1.0149650831869401, 1.0149650831869401, 0.},
     {0., -0.09544592783001707, 0.},
     {1.2508164096106962, 0.3338069521032678, -3.0649628355988305}}

    For Python/Numpy:

    >>> a, A = RandomPH(3, 8, 4)
    >>> print(a)
    [[ 0.50551  0.       0.49449]]
    >>> print(A)
    [[-2.40914  0.91736  0.74455]
     [ 1.33012 -2.52955  1.19943]
     [ 0.       0.85372 -0.85372]]

