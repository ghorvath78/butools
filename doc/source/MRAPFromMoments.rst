butools.map.MRAPFromMoments
===========================

.. currentmodule:: butools.map

.. np:function:: MRAPFromMoments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`H = MRAPFromMoments(moms, Nm)`
        * - Mathematica:
          - :code:`H = MRAPFromMoments[moms, Nm]`
        * - Python/Numpy:
          - :code:`H = MRAPFromMoments(moms, Nm)`

    Creates a marked rational arrival process that has the same 
    marginal and lag-1 joint moments as given (see [1]_).

    Parameters
    ----------
    moms : vector of doubles
        The list of marginal moments. To obtain a marked 
        rational process of order M, 2*M-1 marginal moments
        are required.
    Nm : list of matrices, shape (M,M)
        The list of lag-1 joint moment matrices. The 
        length of the list determines K, the number of arrival 
        types of the rational process.
    
    Returns
    -------
    H : list of matrices, shape (M,M)
        The H0, H1, ..., HK matrices of the marked rational 
        process
    
    Notes
    -----
    There is no guarantee that the returned matrices define
    a valid stochastic process. The joint densities may be
    negative.

    References
    ----------
    .. [1] Andras Horvath, Gabor Horvath, Miklos Telek, "A 
           traffic based decomposition of two-class queueing
           networks with priority service," Computer Networks 
           53:(8) pp. 1235-1248. (2009)

    Examples
    --------
    For Matlab:
    
    >>> moms = [0.99805, 2.0262, 6.2348, 25.738, 133.3];
    >>> jmoms1 = [0.34, 0.34055, 0.69322; 0.34966, 0.35186, 0.71829; 0.72396, 0.73079, 1.4947];
    >>> jmoms2 = [0.29553, 0.27734, 0.54078; 0.29028, 0.27211, 0.53013; 0.58391, 0.54682, 1.0646];
    >>> jmoms3 = [0.36447, 0.38016, 0.79224; 0.35811, 0.37446, 0.78148; 0.71836, 0.75239, 1.5716];
    >>> H = MRAPFromMoments(moms,{jmoms1,jmoms2,jmoms3});
    >>> H{1}
          -1.9473       3.0344      -2.1704
         -0.33434     -0.88118      0.21355
         -0.33363      0.21321     -0.88152
    >>> H{2}
          0.22304       7.4962      -7.3973
         0.092805       7.9131      -7.6658
         0.086985       7.8804      -7.6261
    >>> H{3}
          0.49395      0.92951      -1.0798
           0.4418      -2.4454       2.2995
          0.44568      -2.3908       2.2415
    >>> H{4}
         -0.14271      -5.4356       5.9959
         -0.19985      -5.1102       5.6761
         -0.19945      -5.1791       5.7429
    >>> MarginalMomentsFromMRAP(H,5,1e-11)
          0.99805       2.0262       6.2348       25.738        133.3
    >>> rjmoms = LagkJointMomentsFromMRAP(H,2,1,1e-11);
    >>> rjmoms{1}
             0.34      0.34055      0.69322
          0.34966      0.35186      0.71829
          0.72396      0.73079       1.4947
    >>> rjmoms{2}
          0.29553      0.27734      0.54078
          0.29028      0.27211      0.53013
          0.58391      0.54682       1.0646
    >>> rjmoms{3}
          0.36447      0.38016      0.79224
          0.35811      0.37446      0.78148
          0.71836      0.75239       1.5716
    
    For Python/Numpy:
    
    >>> moms = [0.99804861334547901, 2.0262345014785748, 6.2347773609141584, 25.738432370871866, 133.30214790750446]
    >>> jmoms1 = ml.matrix([[ 0.34000092,  0.34054528,  0.69321737],[ 0.34966003,  0.35185598,  0.71829122],[ 0.72396246,  0.73078643,  1.49469082]])
    >>> jmoms2 = ml.matrix([[ 0.2955327 ,  0.27734236,  0.54077735],[ 0.29028126,  0.27210886,  0.53013258],[ 0.5839079 ,  0.546818  ,  1.06457498]])
    >>> jmoms3 = ml.matrix([[ 0.36446638,  0.38016098,  0.79223978],[ 0.35810733,  0.37445996,  0.78147648],[ 0.71836414,  0.75238761,  1.57164577]])
    >>> H = MRAPFromMoments(moms,(jmoms1,jmoms2,jmoms3))
    >>> print(H[0])
    [[-1.94730519  3.03441393 -2.17035543]
     [-0.3343442  -0.88117688  0.21354779]
     [-0.33362593  0.21320675 -0.88151792]]
    >>> print(H[1])
    [[ 0.22303552  7.49252265 -7.3936431 ]
     [ 0.09280522  7.91365871 -7.66639786]
     [ 0.08698437  7.87984763 -7.6255667 ]]
    >>> print(H[2])
    [[ 0.49394577  0.92090262 -1.071202  ]
     [ 0.44180109 -2.44500194  2.29906455]
     [ 0.44567871 -2.39209787  2.24277651]]
    >>> print(MarginalMomentsFromMRAP(H,5,1e-7))
    [0.99804866709421602, 2.0262346537338147, 6.2347778420947577, 25.738434141923204, 133.30215529953227]
    >>> rjmoms = LagkJointMomentsFromMRAP(H,2,1,1e-7)
    >>>  print(rjmoms[0])
    [[ 0.34000159  0.34054613  0.69321931]
     [ 0.34966065  0.35185676  0.71829299]
     [ 0.72396359  0.73078784  1.49469405]]
     >>> print(rjmoms[1])
    [[ 0.29553296  0.27734252  0.54077755]
     [ 0.29028152  0.27210903  0.5301328 ]
     [ 0.58390841  0.54681833  1.06457543]]
    >>> print(rjmoms[2])
    [[ 0.36446544  0.38016002  0.7922378 ]
     [ 0.3581065   0.37445911  0.78147472]
     [ 0.71836265  0.75238608  1.57164262]]





    

