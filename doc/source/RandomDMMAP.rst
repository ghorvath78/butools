butools.dmap.RandomDMMAP
========================

.. currentmodule:: butools.dmap

.. np:function:: RandomDMMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`D = RandomDMMAP(order, types, mean, zeroEntries, maxTrials, prec)`
        * - Mathematica:
          - :code:`D = RandomDMMAP[order, types, mean, zeroEntries, maxTrials, prec]`
        * - Python/Numpy:
          - :code:`D = RandomDMMAP(order, types, mean, zeroEntries, maxTrials, prec)`

    Returns a random discrete Markovian arrival process.

    Parameters
    ----------
    order : int
        The size of the DMAP
    mean : double, optional
        The mean inter-arrival times of the DMMAP
    types : int
        The number of different arrival types
    zeroEntries : int, optional
        The number of zero entries in the D0 and D1 matrices
    maxTrials : int, optional
        The maximum number of trials to find a proper DMMAP 
        (that has an irreducible phase process and none of 
        its parameters is all-zero)
    prec : double, optional
        Numerical precision for checking the irreducibility.
        The default value is 1e-14.

    Returns
    -------
    D : list/cell of matrices of shape(M,M), length(types+1)
        The D0...Dtypes matrices of the DMMAP 

    Notes
    -----
    If it fails, try to increase the 'maxTrials' parameter,
    or/and the 'mean' parameter.

    Examples
    --------
    For Matlab:
    
    >>> D=RandomDMMAP(4,2,5.62,10);
    >>> D{1}
          0.75138     0.017781    0.0064893      0.03438
         0.018164      0.70408    0.0051205     0.052919
         0.012098            0       0.7573            0
                0     0.007301            0        0.863
    >>> D{2}
         0.012518     0.027781      0.02569     0.031729
         0.039446            0            0     0.073397
         0.037221    0.0013345     0.056364      0.04962
          0.02829            0     0.050616            0
    >>> D{3}
         0.031034     0.032878     0.028048   0.00029258
         0.035816     0.058629            0      0.01243
                0      0.03777    0.0081643     0.040125
         0.020318     0.025133    0.0014527    0.0038849
    >>> mean = MarginalMomentsFromDMMAP(D,1)
             5.62
    >>> sum(sum((D{1}==0) + (D{2}==0) + (D{3}==0)))
        10                 

    For Python/Numpy:
    
    >>> D=RandomDMMAP(4,3,5.62,10)
    >>> print(D[0])
    [[ 0.75291635  0.02838976  0.00623525  0.00355369]
     [ 0.01363156  0.88103764  0.00648297  0.00752309]
     [ 0.05029372  0.04629122  0.5646528   0.03090788]
     [ 0.03333324  0.12876645  0.          0.42277574]]
    >>> print(D[1])
    [[ 0.0229293   0.          0.03024724  0.01160849]
     [ 0.00811442  0.0006299   0.01401004  0.0117653 ]
     [ 0.00889419  0.          0.01204694  0.        ]
     [ 0.          0.10835255  0.05132838  0.03689629]]
    >>> print(D[2])
    [[ 0.00836467  0.0221283   0.01718271  0.02484342]
     [ 0.00497048  0.00130812  0.00556127  0.01410468]
     [ 0.03284471  0.04024043  0.03279771  0.023646  ]
     [ 0.          0.          0.08162577  0.01997635]]
    >>> print(D[3])
    [[ 0.02548598  0.01142023  0.01925548  0.01543914]
     [ 0.          0.0133691   0.00716768  0.01032375]
     [ 0.04677512  0.04188656  0.03148354  0.0372392 ]
     [ 0.          0.03594718  0.08099805  0.        ]]
    >>> m = MarginalMomentsFromDMMAP(D,1)[0]
    >>> print(m)
    5.62
    >>> print(np.sum(D[0]==0)+np.sum(D[1]==0)+np.sum(D[2]==0)+np.sum(D[3]==0))
    10
    
