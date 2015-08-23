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
    ========
    For Matlab:

    >>> G0 = [-1.05, 0.03, 0.07; 0.19, -1.63, 0.06; 0, 0.2, -1.03];
    >>> G1 = [0.16, 0.11, 0; 0.1, 0.16, 0; 0.27, 0, 0.19];
    >>> G2 = [0.01, 0.09, 0.13; 0.26, 0.21, 0.05; 0, 0.16, 0.07];
    >>> G3 = [0.19, 0.06, 0.2; 0.17, 0.16, 0.27; 0, 0, 0.14];
    >>> G = {G0, G1, G2, G3};
    >>> moms = MarginalMomentsFromMRAP(G, 5);
    >>> disp(moms);
          0.99805       2.0262       6.2348       25.738        133.3
    >>> Nm = LagkJointMomentsFromMRAP(G, 2, 1);
    >>> [Nm1, Nm2, Nm3] = Nm{:};
    >>> disp(Nm1);
             0.34      0.34055      0.69322
          0.34966      0.35186      0.71829
          0.72396      0.73079       1.4947
    >>> disp(Nm2);
          0.29553      0.27734      0.54078
          0.29028      0.27211      0.53013
          0.58391      0.54682       1.0646
    >>> disp(Nm3);
          0.36447      0.38016      0.79224
          0.35811      0.37446      0.78148
          0.71836      0.75239       1.5716
    >>> H = MRAPFromMoments(moms, Nm);
    >>> disp(H{1});
          -1.9473       3.0344      -2.1704
         -0.33434     -0.88118      0.21355
         -0.33363      0.21321     -0.88152
    >>> disp(H{2});
          0.22304       7.4962      -7.3973
         0.092805       7.9131      -7.6658
         0.086985       7.8804      -7.6261
    >>> disp(H{3});
          0.49395      0.92951      -1.0798
           0.4418      -2.4454       2.2995
          0.44568      -2.3908       2.2415
    >>> disp(H{4});
         -0.14271      -5.4356       5.9959
         -0.19985      -5.1102       5.6761
         -0.19945      -5.1791       5.7429
    >>> BuToolsCheckPrecision = 10.^-10;
    >>> rmoms = MarginalMomentsFromMRAP(H, 5);
    >>> disp(rmoms);
          0.99805       2.0262       6.2348       25.738        133.3
    >>> rNm = LagkJointMomentsFromMRAP(H, 2, 1);
    >>> [rNm1, rNm2, rNm3] = rNm{:};
    >>> disp(rNm1);
             0.34      0.34055      0.69322
          0.34966      0.35186      0.71829
          0.72396      0.73079       1.4947
    >>> disp(rNm2);
          0.29553      0.27734      0.54078
          0.29028      0.27211      0.53013
          0.58391      0.54682       1.0646
    >>> disp(rNm3);
          0.36447      0.38016      0.79224
          0.35811      0.37446      0.78148
          0.71836      0.75239       1.5716

    For Mathematica:

    
    For Python/Numpy:

    >>> G0 = ml.matrix([[-1.05, 0.03, 0.07],[0.19, -1.63, 0.06],[0, 0.2, -1.03]])
    >>> G1 = ml.matrix([[0.16, 0.11, 0],[0.1, 0.16, 0],[0.27, 0, 0.19]])
    >>> G2 = ml.matrix([[0.01, 0.09, 0.13],[0.26, 0.21, 0.05],[0, 0.16, 0.07]])
    >>> G3 = ml.matrix([[0.19, 0.06, 0.2],[0.17, 0.16, 0.27],[0, 0, 0.14]])
    >>> G = [G0, G1, G2, G3]
    >>> moms = MarginalMomentsFromMRAP(G, 5)
    >>> print(moms)
    [0.99804861334547901, 2.0262345014785748, 6.2347773609141584, 25.738432370871866, 133.30214790750446]
    >>> Nm = LagkJointMomentsFromMRAP(G, 2, 1)
    >>> Nm1, Nm2, Nm3 = Nm
    >>> print(Nm1)
    [[ 0.34     0.34055  0.69322]
     [ 0.34966  0.35186  0.71829]
     [ 0.72396  0.73079  1.49469]]
    >>> print(Nm2)
    [[ 0.29553  0.27734  0.54078]
     [ 0.29028  0.27211  0.53013]
     [ 0.58391  0.54682  1.06457]]
    >>> print(Nm3)
    [[ 0.36447  0.38016  0.79224]
     [ 0.35811  0.37446  0.78148]
     [ 0.71836  0.75239  1.57165]]
    >>> H = MRAPFromMoments(moms, Nm)
    >>> print(H[0])
    [[-1.94731  3.03441 -2.17036]
     [-0.33434 -0.88118  0.21355]
     [-0.33363  0.21321 -0.88152]]
    >>> print(H[1])
    [[ 0.22304  7.49615 -7.39728]
     [ 0.0928   7.91308 -7.66582]
     [ 0.08698  7.8804  -7.62612]]
    >>> print(H[2])
    [[ 0.49395  0.92951 -1.07981]
     [ 0.4418  -2.4454   2.29947]
     [ 0.44568 -2.39078  2.24146]]
    >>> print(H[3])
    [[-0.14271 -5.43556  5.99595]
     [-0.19985 -5.1102   5.67609]
     [-0.19945 -5.17914  5.7429 ]]
    >>> butools.checkPrecision = 10.**-10
    >>> rmoms = MarginalMomentsFromMRAP(H, 5)
    >>> print(rmoms)
    [0.9980486133458466, 2.026234501479784, 6.2347773609187414, 25.738432370892887, 133.30214790761983]
    >>> rNm = LagkJointMomentsFromMRAP(H, 2, 1)
    >>> rNm1, rNm2, rNm3 = rNm
    >>> print(rNm1)
    [[ 0.34     0.34055  0.69322]
     [ 0.34966  0.35186  0.71829]
     [ 0.72396  0.73079  1.49469]]
    >>> print(rNm2)
    [[ 0.29553  0.27734  0.54078]
     [ 0.29028  0.27211  0.53013]
     [ 0.58391  0.54682  1.06457]]
    >>> print(rNm3)
    [[ 0.36447  0.38016  0.79224]
     [ 0.35811  0.37446  0.78148]
     [ 0.71836  0.75239  1.57165]]

