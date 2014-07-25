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
    
    

