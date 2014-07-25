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

