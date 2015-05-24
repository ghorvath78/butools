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
    --------
    For Matlab:
    
    >>> [D0,D1]=RandomDMAP(4,5.62,10);
    >>> D0
           0.6394      0.11806            0     0.022015
         0.075135      0.50301     0.046918      0.15842
                0     0.080273      0.73782     0.056975
                0      0.19654      0.17778      0.36056
    >>> D1
         0.067826            0     0.075489      0.07721
         0.088654            0            0      0.12786
                0    0.0049056     0.074522     0.045509
                0            0      0.26511            0
    >>> mean = MarginalMomentsFromDMAP(D0,D1,1)
             5.62
    >>> sum(sum((D0==0) + (D1==0)))
        10                 

    For Python/Numpy:
    
    >>> [D0,D1]=RandomDMAP(4,5.62,10)
    >>> print(D0)
    [[ 0.88058332  0.00922074  0.          0.01169603]
     [ 0.17186087  0.13291209  0.151581    0.08091293]
     [ 0.          0.09170655  0.36737607  0.05133128]
     [ 0.12319488  0.          0.          0.33909318]]
    >>> print(D1)
    [[ 0.02226474  0.03165048  0.01921701  0.02536769]
     [ 0.11526746  0.          0.18367335  0.1637923 ]
     [ 0.23220325  0.          0.25738285  0.        ]
     [ 0.53771194  0.          0.          0.        ]]
    >>> m = MarginalMomentsFromDMAP(D0,D1,1)[0]
    >>> print(m)
    5.62
    >>> print(np.sum((D0==0))+np.sum((D1==0)))
    10
    
 
    
