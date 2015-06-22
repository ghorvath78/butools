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
         1     0     0
    >>> disp(A);
         -0.28487      0.11715    0.0062532
         0.031733     -0.16649     0.017921
         0.076758            0    -0.076758

    For Mathematica:

    >>> {a,A} = RandomPH[3,8,4];
    >>> Print[a];
    {0, 0, 1}
    >>> Print[A];
    {{-0.6999227429898282, 0.3967535668149389, 0.028566129628349583},
     {0.446826923286755, -0.9416070080107307, 0.38815288975615425},
     {0.45083462621244397, 0., -0.45083462621244397}}

    For Python/Numpy:

    >>> a,A = RandomPH(3,8,4)
    >>> print(a)
    [[ 0.       0.47702  0.52298]]
    >>> print(A)
    [[-4.64155  1.23652  2.35056]
     [ 2.35868 -3.01819  0.6595 ]
     [ 0.47849  0.      -0.47849]]

