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

    >>> [a,A] = RandomDPH(3,10,5);
    CheckProbMatrix: the matrix has negative element (precision: 1e-07)!
    >>> disp(a);
          0.35022      0.64978            0
    >>> disp(A);
          0.59299            0      0.14264
          0.31704      0.45572      0.22724
                0      0.56445      0.43555

    For Mathematica:

    >>> {a,A} = RandomDPH[3,10,5];
    "CheckProbMatrix: the matrix has negative element (at precision "1.*^-7")!"
    >>> Print[a];
    {1, 0, 0}
    >>> Print[A];
    {{0.6747139991829699, 0., 0.06813341255217867},
     {0.31858283653858926, 0.42042758912527506, 0.03433200607772197},
     {0.021164290306528725, 0.021429694043095273, 0.957406015650376}}

    For Python/Numpy:

    >>> a,A = RandomDPH(3,10,5)
    >>> print(a)
    [[ 1.  0.  0.]]
    >>> print(A)
    [[ 0.84936  0.07094  0.00315]
     [ 0.0107   0.82367  0.     ]
     [ 0.10697  0.       0.89303]]

