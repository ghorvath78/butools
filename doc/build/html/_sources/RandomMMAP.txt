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
    ========
    For Matlab:

    >>> Dm = RandomMMAP(4, 3, 1.62, 10);
    >>> disp(Dm{1});
         -0.84147     0.066601     0.054766            0
         0.071719     -0.67551            0     0.016834
         0.016327      0.11655     -0.57318     0.043691
         0.095968     0.079712            0     -0.87041
    >>> disp(Dm{2});
         0.085267      0.12324     0.018852      0.08092
          0.04971     0.029681     0.015256      0.10755
          0.04997            0            0    0.0019056
         0.098333     0.062911      0.06319     0.055274
    >>> disp(Dm{3});
          0.11299     0.092558            0     0.028328
         0.043418      0.10553    0.0093983     0.063497
          0.10159     0.071162            0    0.0073957
                0     0.053923      0.04684     0.060224
    >>> disp(Dm{4});
          0.10217     0.055739    0.0096715      0.01037
                0     0.017933     0.076958     0.068026
         0.011933    0.0053227     0.029046      0.11829
           0.0379      0.10115            0      0.11498
    >>> m = MarginalMomentsFromMMAP(Dm, 1);
    >>> disp(m);
             1.62

    For Mathematica:

    
    For Python/Numpy:

    >>> Dm = RandomMMAP(4, 3, 1.62, 10)
    >>> print(Dm[0])
    [[-0.77004  0.07746  0.00752  0.12555]
     [ 0.00951 -0.96329  0.11107  0.04652]
     [ 0.       0.07548 -0.66821  0.03229]
     [ 0.12584  0.       0.03691 -0.77768]]
    >>> print(Dm[1])
    [[ 0.05194  0.03399  0.00531  0.04182]
     [ 0.02984  0.07001  0.1335   0.06569]
     [ 0.03889  0.04339  0.       0.     ]
     [ 0.06669  0.       0.10119  0.0133 ]]
    >>> print(Dm[2])
    [[ 0.08729  0.03552  0.03436  0.0247 ]
     [ 0.11078  0.02114  0.       0.02779]
     [ 0.11427  0.06711  0.0964   0.12043]
     [ 0.05964  0.1119   0.       0.06175]]
    >>> print(Dm[3])
    [[ 0.07868  0.04748  0.02386  0.09457]
     [ 0.14246  0.       0.07483  0.12014]
     [ 0.       0.05206  0.00753  0.02037]
     [ 0.0462   0.04298  0.11128  0.     ]]
    >>> m = MarginalMomentsFromMMAP(Dm, 1)[0]
    >>> print(m)
    1.62

