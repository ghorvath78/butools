butools.dmap.DRAPFromMoments
============================

.. currentmodule:: butools.dmap

.. np:function:: DRAPFromMoments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[H0, H1] = DRAPFromMoments(moms, Nm)`
        * - Mathematica:
          - :code:`{H0, H1} = DRAPFromMoments[moms, Nm]`
        * - Python/Numpy:
          - :code:`H0, H1 = DRAPFromMoments(moms, Nm)`

    Creates a discrete rational arrival process that has the 
    same marginal and lag-1 joint moments as given (see [1]_).

    Parameters
    ----------
    moms : vector of doubles
        The list of marginal moments. To obtain a rational 
        process of order M, 2*M-1 marginal moments are 
        required.
    Nm : matrix, shape (M,M)
        The matrix of lag-1 joint moments. 
    
    Returns
    -------
    H0 : matrix, shape (M,M)
        The H0 matrix of the discrete rational process
    H1 : matrix, shape (M,M)
        The H1 matrix of the discrete rational process
    
    References
    ----------
    .. [1] G Horvath, M Telek, "A minimal representation of 
           Markov arrival processes and a moments matching 
           method," Performance Evaluation 64:(9-12) pp. 
           1153-1168. (2007)       

    Examples
    ========
    For Matlab:

    >>> G0 = [0, 0.02, 0, 0; 0, 0.17, 0.2, 0.14; 0.16, 0.17, 0.02, 0.18; 0, 0, 0, 0.12];
    >>> G1 = [0, 0.88, 0.1, 0; 0.18, 0.07, 0.14, 0.1; 0.13, 0.15, 0.15, 0.04; 0.31, 0.18, 0.12, 0.27];
    >>> moms = MarginalMomentsFromDRAP(G0, G1, 5);
    >>> disp(moms);
           1.4955       2.9542       7.8852       27.282       116.17
    >>> Nm = LagkJointMomentsFromDRAP(G0, G1, 2, 1);
    >>> disp(Nm);
                1       1.4955       2.9542
           1.4955       2.2037       4.2827
           2.9542       4.2875       8.1899
    >>> [H0, H1] = DRAPFromMoments(moms, Nm);
    >>> disp(H0);
          0.56447      0.47188     -0.69474
         -0.50857     -0.10551      0.95921
          0.18477      0.26121     -0.13431
    >>> disp(H1);
           2.3994       1.1243      -2.8653
          -1.7535     -0.59009       2.9984
          0.95074      0.51879      -0.7812
    >>> rmoms = MarginalMomentsFromDRAP(H0, H1, 5);
    >>> disp(rmoms);
           1.4955       2.9542       7.8852       27.282       116.17
    >>> rNm = LagkJointMomentsFromDRAP(H0, H1, 2, 1);
    >>> disp(rNm);
                1       1.4955       2.9542
           1.4955       2.2037       4.2827
           2.9542       4.2875       8.1899

    For Mathematica:

    >>> G0 = {{0, 0.02, 0, 0},{0, 0.17, 0.2, 0.14},{0.16, 0.17, 0.02, 0.18},{0, 0, 0, 0.12}};
    >>> G1 = {{0, 0.88, 0.1, 0},{0.18, 0.07, 0.14, 0.1},{0.13, 0.15, 0.15, 0.04},{0.31, 0.18, 0.12, 0.27}};
    >>> moms = MarginalMomentsFromDRAP[G0, G1, 5];
    >>> Print[moms];
    {1.4955358592094412, 2.9542479654368474, 7.885226907678561, 27.282328108669493, 116.17171481905851}
    >>> Nm = LagkJointMomentsFromDRAP[G0, G1, 2, 1];
    >>> Print[Nm];
    {{1., 1.4955358592094412, 2.954247965436847},
     {1.4955358592094414, 2.2037182406034797, 4.282673397390962},
     {2.9542479654368474, 4.287487747878976, 8.189899409259828}}
    >>> {H0, H1} = DRAPFromMoments[moms, Nm];
    >>> Print[H0];
    {{0.5644738962225417, 0.47187846354848406, -0.6947446288880126},
     {-0.5085686970022437, -0.10550993297233946, 0.9592122034338286},
     {0.18477321176067846, 0.261205874785389, -0.13431385414476443}}
    >>> Print[H1];
    {{2.3993787625484853, 1.1243091633982, -2.865295656829698},
     {-1.7534589569688064, -0.5900943661745801, 2.9984197496841354},
     {0.9507424301490942, 0.5187877102745091, -0.7811953728249019}}
    >>> rmoms = MarginalMomentsFromDRAP[H0, H1, 5];
    >>> Print[rmoms];
    {1.495535859209443, 2.954247965436855, 7.885226907678592, 27.28232810866962, 116.17171481905912}
    >>> rNm = LagkJointMomentsFromDRAP[H0, H1, 2, 1];
    >>> Print[rNm];
    {{0.9999999999999997, 1.4955358592094425, 2.954247965436854},
     {1.4955358592094437, 2.2037182406034845, 4.282673397390974},
     {2.9542479654368594, 4.287487747878995, 8.189899409259864}}

    For Python/Numpy:

    >>> G0 = ml.matrix([[0, 0.02, 0, 0],[0, 0.17, 0.2, 0.14],[0.16, 0.17, 0.02, 0.18],[0, 0, 0, 0.12]])
    >>> G1 = ml.matrix([[0, 0.88, 0.1, 0],[0.18, 0.07, 0.14, 0.1],[0.13, 0.15, 0.15, 0.04],[0.31, 0.18, 0.12, 0.27]])
    >>> moms = MarginalMomentsFromDRAP(G0, G1, 5)
    >>> print(moms)
    [1.4955358592094412, 2.9542479654368474, 7.885226907678561, 27.282328108669493, 116.17171481905851]
    >>> Nm = LagkJointMomentsFromDRAP(G0, G1, 2, 1)
    >>> print(Nm)
    [[ 1.       1.49554  2.95425]
     [ 1.49554  2.20372  4.28267]
     [ 2.95425  4.28749  8.1899 ]]
    >>> H0, H1 = DRAPFromMoments(moms, Nm)
    >>> print(H0)
    [[ 0.56447  0.47188 -0.69474]
     [-0.50857 -0.10551  0.95921]
     [ 0.18477  0.26121 -0.13431]]
    >>> print(H1)
    [[ 2.39938  1.12431 -2.8653 ]
     [-1.75346 -0.59009  2.99842]
     [ 0.95074  0.51879 -0.7812 ]]
    >>> rmoms = MarginalMomentsFromDRAP(H0, H1, 5)
    >>> print(rmoms)
    [1.495535859209453, 2.9542479654368994, 7.885226907678768, 27.282328108670363, 116.17171481906257]
    >>> rNm = LagkJointMomentsFromDRAP(H0, H1, 2, 1)
    >>> print(rNm)
    [[ 1.       1.49554  2.95425]
     [ 1.49554  2.20372  4.28267]
     [ 2.95425  4.28749  8.1899 ]]

