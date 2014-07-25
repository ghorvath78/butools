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
    --------
    For Matlab:
    
    >>> [D0,D1]=RandomMAP(4,1.62,10);
    >>> D0
         -0.87934            0            0            0
                0      -1.2074      0.42836      0.25102
          0.11181            0     -0.69366      0.24076
          0.44849      0.42844      0.38088       -2.319
    >>> D1
          0.38098      0.22325            0      0.27511
                0      0.39026            0      0.13776
                0            0      0.18773      0.15335
          0.12257      0.21414      0.53865      0.18586
    >>> mean = MarginalMomentsFromMAP(D0,D1,1)
             1.62
    >>> sum(sum((D0==0) + (D1==0)))
        10                 

