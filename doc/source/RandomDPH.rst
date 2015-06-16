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
    >>> disp(a);
         1     0     0
    >>> disp(A);
          0.33352      0.29938      0.11926
                0      0.93251     0.042746
                0            0      0.74869

    For Mathematica:

    >>> {a,A} = RandomDPH[3,10,5];
    >>> Print[a];
    {0, 1, 0}
    >>> Print[A];
    {{0.634314117595034, 0.03145872350022707, 0.163933538385891},
     {0.22501882446572802, 0.7742605331911611, 0.},
     {0., 0.0740389003080393, 0.6085494083431278}}

    For Python/Numpy:

    >>> a,A = RandomDPH(3,10,5)
    >>> print(a)
    [[ 1.  0.  0.]]
    >>> print(A)
    [[ 0.80329  0.05712  0.05   ]
     [ 0.10619  0.84325  0.     ]
     [ 0.       0.       0.76827]]

