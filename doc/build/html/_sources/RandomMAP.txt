butools.map.RandomMAP
=====================

.. currentmodule:: butools.map

.. np:function:: RandomMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[D0, D1] = RandomMAP(order, mean, zeroEntries, maxTrials, prec)`
        * - Mathematica:
          - :code:`{D0, D1} = RandomMAP[order, mean, zeroEntries, maxTrials, prec]`
        * - Python/Numpy:
          - :code:`D0, D1 = RandomMAP(order, mean, zeroEntries, maxTrials, prec)`

    Returns a random Markovian arrival process with given mean 
    value.

    Parameters
    ----------
    order : int
        The size of the MAP
    mean : double, optional
        The mean inter-arrival times of the MAP
    zeroEntries : int, optional
        The number of zero entries in the D0 and D1 matrices
    maxTrials : int, optional
        The maximum number of trials to find a proper MAP 
        (that has an irreducible phase process and none of 
        its parameters is all-zero)
    prec : double, optional
        Numerical precision for checking the irreducibility.
        The default value is 1e-14.

    Returns
    -------
    D0 : vector, shape (1,M)
        The D0 matrix of the MAP
    D1 : matrix, shape (M,M)
        The D1 matrix of the MAP

    Examples
    ========
    For Matlab:

    >>> [D0, D1] = RandomMAP(4, 1.62, 10);
    >>> disp(D0);
         -0.93599            0     0.037316            0
                0      -2.1684      0.36013      0.64269
          0.77567            0      -2.0228            0
          0.44227            0      0.11204     -0.67496
    >>> disp(D1);
         0.078503     0.025728      0.22381      0.57063
                0       0.6187      0.39583      0.15103
          0.40278      0.47302       0.2751     0.096175
          0.12066            0            0            0
    >>> m = MarginalMomentsFromMAP(D0, D1, 1);
    >>> disp(m);
             1.62

    For Mathematica:

    >>> {D0, D1} = RandomMAP[4, 1.62, 10];
    >>> Print[D0];
    {{-1.2603255570025995, 0., 1.2603255570025995, 0.},
     {1.0421148823453172, -3.0725926090541065, 0., 1.2512140562551135},
     {0.0729686946916973, 1.1513050282525865, -3.021511605493252, 0.8603731650240692},
     {0.905645182073805, 0., 0.3797626274429608, -2.652481770444762}}
    >>> Print[D1];
    {{0., 0., 0., 0.},
     {0.0035497683267016345, 0.025589968292695654, 0.1383759482229132, 0.6117479856113653},
     {0.4204286154800259, 0.36476933475507295, 0.14415600056384018, 0.007510766725959241},
     {0.1266680865059269, 1.2404058744220698, 0., 0.}}
    >>> m = MarginalMomentsFromMAP[D0, D1, 1][[1]];
    >>> Print[m];
    1.6199999999999997

    For Python/Numpy:

    >>> D0, D1 = RandomMAP(4, 1.62, 10)
    >>> print(D0)
    [[-2.72432  0.64111  0.42238  0.     ]
     [ 0.      -1.75493  0.18555  0.48421]
     [ 0.62751  0.      -2.14775  0.48912]
     [ 0.       0.       0.01388 -0.3889 ]]
    >>> print(D1)
    [[ 0.19832  0.41384  0.46643  0.58223]
     [ 0.45159  0.15916  0.       0.47442]
     [ 0.43032  0.02302  0.       0.57777]
     [ 0.       0.37502  0.       0.     ]]
    >>> m = MarginalMomentsFromMAP(D0, D1, 1)[0]
    >>> print(m)
    1.62

