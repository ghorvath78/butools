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
          0.58888      0.41112            0
    >>> disp(A);
           0.9048      0.04125     0.031817
                0       0.6822       0.3178
                0            0      0.49751

    For Mathematica:

    >>> {a,A} = RandomDPH[3,10,5];
    "CheckProbMatrix: the matrix has negative element (at precision "1.*^-7")!"
    "CheckProbMatrix: the matrix has negative element (at precision "1.*^-7")!"
    "CheckProbMatrix: the matrix has negative element (at precision "1.*^-7")!"
    "CheckProbMatrix: the matrix has negative element (at precision "1.*^-7")!"
    "CheckProbMatrix: the matrix has negative element (at precision "1.*^-7")!"
    >>> Print[a];
    {0, 1, 0}
    >>> Print[A];
    {{0.9978146840537366, 0., 0.0018115430265152948},
     {0.0002596368902400116, 0.8783862382976336, 0.005417700244985481},
     {0., 0.07246410696945481, 0.9001849760112214}}

    For Python/Numpy:

    >>> a,A = RandomDPH(3,10,5)
    CheckProbMatrix: the matrix has negative element (precision: 1e-07)!
    >>> print(a)
    [[ 1.  0.  0.]]
    >>> print(A)
    [[ 0.85043  0.00248  0.13613]
     [ 0.01293  0.87899  0.     ]
     [ 0.       0.       0.71197]]

