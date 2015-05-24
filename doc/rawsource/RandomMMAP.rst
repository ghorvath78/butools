butools.map.RandomMMAP
======================

.. currentmodule:: butools.map

.. np:function:: RandomMMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`D = RandomMMAP(order, types, mean, zeroEntries, maxTrials, prec)`
        * - Mathematica:
          - :code:`D = RandomMMAP[order, types, mean, zeroEntries, maxTrials, prec]`
        * - Python/Numpy:
          - :code:`D = RandomMMAP(order, types, mean, zeroEntries, maxTrials, prec)`

    Returns a random Markovian arrival process with given mean 
    value.

    Parameters
    ----------
    order : int
        The size of the MAP
    types : int
        The number of different arrival types
    mean : double, optional
        The mean inter-arrival times of the MMAP
    zeroEntries : int, optional
        The number of zero entries in the D0 and D1 matrices
    maxTrials : int, optional
        The maximum number of trials to find a proper MMAP 
        (that has an irreducible phase process and none of 
        its parameters is all-zero)
    prec : double, optional
        Numerical precision for checking the irreducibility.
        The default value is 1e-14.

    Returns
    -------
    D : list/cell of matrices of shape(M,M), length(types+1)
        The D0...Dtypes matrices of the MMAP 

    Examples
    --------
    For Matlab:
    
    >>> D=RandomMMAP(4,2,1.62,10);
    >>> D{1}
          -1.0886      0.12136      0.22192      0.16842
         0.014734     -0.95569      0.12078       0.1985
          0.04786            0     -0.81391      0.03485
                0            0      0.19784     -0.62705
    >>> D{2}
                0      0.21883      0.11589    0.0065551
          0.10771     0.046281            0     0.026041
         0.018566     0.077934      0.11477     0.055241
         0.021011     0.023357      0.13133            0
    >>> D{3}
                0      0.14526     0.090391            0
                0       0.1356       0.2126     0.093445
         0.013967      0.23254      0.15518     0.063006
                0     0.077494      0.16397     0.012044
    >>> mean = MarginalMomentsFromMMAP(D,1)
             1.62
    >>> sum(sum((D{1}==0) + (D{2}==0) + (D{3}==0)))
        10                 

    For Python/Numpy:

    >>> D=RandomMMAP(4,2,1.62,10)
    >>> print(D[0])
    [[-1.20043764  0.05686545  0.17894297  0.03642876]
     [ 0.0746768  -1.042511    0.05897844  0.03616775]
     [ 0.10703581  0.15982011 -0.78408293  0.        ]
     [ 0.          0.          0.         -0.40896525]]
    >>> print(D[1])
    [[ 0.06671026  0.09825053  0.20481737  0.18312201]
     [ 0.00370989  0.12024972  0.20318316  0.09044313]
     [ 0.08165806  0.05115381  0.1339248   0.        ]
     [ 0.07774122  0.09463614  0.          0.        ]]
    >>> print(D[2])
    [[ 0.11557699  0.06062982  0.14678981  0.05230366]
     [ 0.06895695  0.06299849  0.15758806  0.16555862]
     [ 0.14566491  0.02460299  0.          0.08022243]
     [ 0.          0.0533731   0.18321478  0.        ]]
    >>> print(MarginalMomentsFromMMAP(D,1))
    [1.6200000000000001]
    >>> print(np.sum(D[0]==0) + np.sum(D[1]==0) + np.sum(D[2]==0))
    10

