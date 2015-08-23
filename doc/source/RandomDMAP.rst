butools.dmap.RandomDMAP
=======================

.. currentmodule:: butools.dmap

.. np:function:: RandomDMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[D0, D1] = RandomDMAP(order, mean, zeroEntries, maxTrials, prec)`
        * - Mathematica:
          - :code:`{D0, D1} = RandomDMAP[order, mean, zeroEntries, maxTrials, prec]`
        * - Python/Numpy:
          - :code:`D0, D1 = RandomDMAP(order, mean, zeroEntries, maxTrials, prec)`

    Returns a random disctere Markovian arrival process.

    Parameters
    ----------
    order : int
        The size of the DMAP
    mean : double, optional
        The mean inter-arrival times of the DMAP
    zeroEntries : int, optional
        The number of zero entries in the D0 and D1 matrices
    maxTrials : int, optional
        The maximum number of trials to find a proper DMAP 
        (that has an irreducible phase process and none of 
        its parameters is all-zero)
    prec : double, optional
        Numerical precision for checking the irreducibility.
        The default value is 1e-14.

    Returns
    -------
    D0 : vector, shape (1,M)
        The D0 matrix of the DMAP
    D1 : matrix, shape (M,M)
        The D1 matrix of the DMAP

    Notes
    -----
    If it fails, try to increase the 'maxTrials' parameter,
    or/and the 'mean' parameter.

    Examples
    ========
    For Matlab:

    >>> [D0, D1] = RandomDMAP(4, 5.62, 10);
    CheckProbMatrix: the matrix has negative element (precision: 1e-12)!
    CheckDMAPRepresentation: D0 isn't a transient probability matrix!
    >>> disp(D0);
          0.32242      0.12926            0       0.1379
          0.11054      0.70995     0.060634     0.012054
                0     0.038601      0.86412     0.021564
                0            0      0.12058      0.42758
    >>> disp(D1);
           0.1503      0.15115    0.0072573       0.1017
         0.069721     0.037105            0            0
         0.035809     0.018712     0.021192            0
                0      0.45185            0            0
    >>> m = MarginalMomentsFromDMAP(D0, D1, 1);
    >>> disp(m);
             5.62

    For Mathematica:

    
    For Python/Numpy:

    >>> D0, D1 = RandomDMAP(4, 5.62, 10)
    CheckDMMAPRepresentation: Some of the matrices D1 ... DM have negative elements!
    >>> print(D0)
    [[ 0.91231  0.       0.       0.02106]
     [ 0.15283  0.16699  0.18905  0.     ]
     [ 0.02497  0.12509  0.52043  0.10254]
     [ 0.0759   0.       0.       0.29181]]
    >>> print(D1)
    [[ 0.00719  0.02013  0.01658  0.02274]
     [ 0.15648  0.       0.18936  0.1453 ]
     [ 0.10408  0.       0.       0.12289]
     [ 0.21349  0.       0.       0.41879]]
    >>> m = MarginalMomentsFromDMAP(D0, D1, 1)[0]
    >>> print(m)
    5.62

