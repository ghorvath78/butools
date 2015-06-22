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

    >>> [D0,D1] = RandomDMAP(4,5.62,10);
    >>> disp(D0);
          0.73286     0.046357     0.034004     0.028796
         0.084185      0.67698            0            0
         0.027256     0.082261      0.74025            0
                0            0      0.80713      0.19287
    >>> disp(D1);
         0.045462     0.020243     0.042927     0.049352
         0.019264      0.09044     0.095789     0.033342
                0      0.10169     0.038362      0.01018
                0            0            0            0
    >>> m = MarginalMomentsFromDMAP(D0,D1,1);
    >>> disp(m);
             5.62

    For Mathematica:

    >>> {D0,D1} = RandomDMAP[4,5.62,10];
    >>> Print[D0];
    {{0.43792221764784234, 0.1315017116128684, 0.015584191529024663, 0.},
     {0., 0.90297335020982, 0., 0.},
     {0.04551398326730008, 0.04581614981871694, 0.3422245239933987, 0.13244936737504273},
     {0.010409491714044225, 0., 0., 0.7212471553509219}}
    >>> Print[D1];
    {{0.16771196633111266, 0.15225054647359362, 0.0033488389139790338, 0.09168052749157933},
     {0., 0.023880459024570835, 0.07314619076560927, 0.},
     {0.12026766413018974, 0.12627847796574468, 0.09656436962438471, 0.09088546382522249},
     {0.06852868252098047, 0.09548585720018578, 0.08366967190544991, 0.020659141308417714}}
    >>> m = MarginalMomentsFromDMAP[D0,D1,1][[1]];
    >>> Print[m];
    5.6200000000000045

    For Python/Numpy:

    >>> D0,D1 = RandomDMAP(4,5.62,10)
    >>> print(D0)
    [[ 0.63757  0.03854  0.05201  0.098  ]
     [ 0.       0.7208   0.05398  0.     ]
     [ 0.       0.02845  0.61792  0.     ]
     [ 0.05187  0.       0.       0.81905]]
    >>> print(D1)
    [[ 0.02668  0.03954  0.10199  0.00567]
     [ 0.08081  0.04892  0.       0.09548]
     [ 0.       0.13362  0.13848  0.08154]
     [ 0.0162   0.       0.       0.11288]]
    >>> m = MarginalMomentsFromDMAP(D0,D1,1)[0]
    >>> print(m)
    5.62

