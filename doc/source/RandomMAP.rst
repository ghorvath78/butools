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

    >>> [D0,D1] = RandomMAP(4,1.62,10);
    >>> disp(D0);
          -2.1564      0.52751      0.23666      0.48227
          0.12958      -1.2882            0      0.40782
         0.095441            0     -0.84224      0.16816
          0.23854            0      0.11448     -0.75312
    >>> disp(D1);
         0.011037      0.54473            0      0.35416
          0.52113      0.22971            0            0
          0.30136            0      0.27727            0
          0.28047      0.11964            0            0
    >>> m = MarginalMomentsFromMAP(D0,D1,1);
    >>> disp(m);
             1.62

    For Mathematica:

    >>> {D0,D1} = RandomMAP[4,1.62,10];
    >>> Print[D0];
    {{-4.951856393334554, 1.4808270612099914, 0.9854373237264185, 0.9601506507976104},
     {0., -0.4862087833155877, 0.4862087833155877, 0.},
     {0., 1.0390511244021063, -3.077603878261726, 0.3851280623523275},
     {1.1400193009034087, 0.9000197303479476, 0.7569716886183007, -7.238543772427738}}
    >>> Print[D1];
    {{1.4661157972439278, 0.023084808049765464, 0., 0.03624075230684058},
     {0., 0., 0., 0.},
     {0., 0., 1.0289936500890735, 0.6244310414182185},
     {1.670741603151747, 0.29891729751978063, 1.3687129632452288, 1.1031611886413246}}
    >>> m = MarginalMomentsFromMAP[D0,D1,1][[1]];
    >>> Print[m];
    1.6199999999999999

    For Python/Numpy:

    >>> D0,D1 = RandomMAP(4,1.62,10)
    >>> print(D0)
    [[-2.12527  0.56416  0.28339  0.13343]
     [ 0.05238 -1.08818  0.       0.45737]
     [ 0.06387  0.      -1.23231  0.     ]
     [ 0.32275  0.       0.58952 -0.91227]]
    >>> print(D1)
    [[ 0.49226  0.30347  0.22972  0.11884]
     [ 0.02567  0.       0.39505  0.15771]
     [ 0.261    0.25525  0.       0.6522 ]
     [ 0.       0.       0.       0.     ]]
    >>> m = MarginalMomentsFromMAP(D0,D1,1)[0]
    >>> print(m)
    1.62

