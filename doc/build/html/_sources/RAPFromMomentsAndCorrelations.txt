butools.map.RAPFromMomentsAndCorrelations
=========================================

.. currentmodule:: butools.map

.. np:function:: RAPFromMomentsAndCorrelations

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[H0, H1] = RAPFromMomentsAndCorrelations(moms, corr)`
        * - Mathematica:
          - :code:`{H0, H1} = RAPFromMomentsAndCorrelations[moms, corr]`
        * - Python/Numpy:
          - :code:`H0, H1 = RAPFromMomentsAndCorrelations(moms, corr)`

    Returns a rational arrival process that has the same moments
    and lag autocorrelation coefficients as given.

    Parameters
    ----------
    moms : vector of doubles
        The vector of marginal moments. To obtain a RAP of 
        size M, 2*M-1 moments are required.
    corr : vector of doubles
        The vector of lag autocorrelation coefficients. To 
        obtain a RAP of size M, 2*M-3 coefficients are needed.
    
    Returns
    -------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational arrival process

    Notes
    -----
    There is no guarantee that the returned matrices define
    a valid stochastic process. The joint densities may be
    negative.

    References
    ----------
    .. [1] Mitchell, Kenneth, and Appie van de Liefvoort. 
           "Approximation models of feed-forward G/G/1/N 
           queueing networks with correlated arrivals." 
           Performance Evaluation 51.2 (2003): 137-152.

    Examples
    ========
    For Matlab:

    >>> H0 = [-6.2, 2., 0; 2., -9., 1.; 1., 0, -3.];
    >>> H1 = [2.2, 0, 2.; 0, 4., 2.; 0, 1., 1.];
    >>> mom = MarginalMomentsFromRAP(H0, H1);
    >>> disp(mom);
          0.29774      0.19284      0.19448      0.26597      0.45833
    >>> corr = LagCorrelationsFromRAP(H0, H1, 3);
    >>> disp(corr);
         0.012394    0.0027412   0.00072384
    >>> [G0, G1] = RAPFromMomentsAndCorrelations(mom, corr);
    >>> disp(G0);
          -8.9629       22.253      -18.544
         -0.99178       -4.667        2.331
          -1.2473       2.4279      -4.5701
    >>> disp(G1);
           2.2027      -1.3173       4.3689
           1.2179       1.8217      0.28809
           1.0212      0.41735        1.951
    >>> rmom = MarginalMomentsFromRAP(G0, G1);
    >>> disp(rmom);
          0.29774      0.19284      0.19448      0.26597      0.45833
    >>> rcorr = LagCorrelationsFromRAP(G0, G1, 3);
    >>> disp(rcorr);
         0.012394    0.0027412   0.00072384

    For Mathematica:

    >>> H0 = {{-6.2, 2., 0},{2., -9., 1.},{1., 0, -3.}};
    >>> H1 = {{2.2, 0, 2.},{0, 4., 2.},{0, 1., 1.}};
    >>> mom = MarginalMomentsFromRAP[H0, H1];
    >>> Print[mom];
    {0.29774127310061604, 0.19283643304803644, 0.19448147792730758, 0.2659732553924553, 0.45833053059627116}
    >>> corr = LagCorrelationsFromRAP[H0, H1, 3];
    >>> Print[corr];
    {0.012393574884970258, 0.0027411959690404088, 0.0007238364213571031}
    >>> {G0, G1} = RAPFromMomentsAndCorrelations[mom, corr];
    >>> Print[G0];
    {{-8.96289388087693, 22.252570107207173, -18.544098091372838},
     {-0.9917815607047362, -4.666992249154709, 2.33103341018933},
     {-1.247298899065379, 2.4279117893845945, -4.570113869959446}}
    >>> Print[G1];
    {{2.2027474563394773, -1.3172514038167056, 4.368925812519816},
     {1.2179262967043782, 1.8217266419760083, 0.2880874609897277},
     {1.0211541630657368, 0.4173538177941065, 1.9509929987803876}}
    >>> rmom = MarginalMomentsFromRAP[G0, G1];
    >>> Print[rmom];
    {0.297741273100616, 0.19283643304803638, 0.19448147792730755, 0.2659732553924553, 0.45833053059627116}
    >>> rcorr = LagCorrelationsFromRAP[G0, G1, 3];
    >>> Print[rcorr];
    {0.012393574884970393, 0.0027411959690408086, 0.0007238364213573696}

    For Python/Numpy:

    >>> H0 = ml.matrix([[-6.2, 2., 0],[2., -9., 1.],[1., 0, -3.]])
    >>> H1 = ml.matrix([[2.2, 0, 2.],[0, 4., 2.],[0, 1., 1.]])
    >>> mom = MarginalMomentsFromRAP(H0, H1)
    >>> print(mom)
    [0.29774127310061604, 0.19283643304803644, 0.19448147792730755, 0.26597325539245531, 0.45833053059627116]
    >>> corr = LagCorrelationsFromRAP(H0, H1, 3)
    >>> print(corr)
    [ 0.01239  0.00274  0.00072]
    >>> G0, G1 = RAPFromMomentsAndCorrelations(mom, corr)
    >>> print(G0)
    [[ -8.96289  22.25257 -18.5441 ]
     [ -0.99178  -4.66699   2.33103]
     [ -1.2473    2.42791  -4.57011]]
    >>> print(G1)
    [[ 2.20275 -1.31725  4.36893]
     [ 1.21793  1.82173  0.28809]
     [ 1.02115  0.41735  1.95099]]
    >>> rmom = MarginalMomentsFromRAP(G0, G1)
    >>> print(rmom)
    [0.29774127310061604, 0.19283643304803638, 0.19448147792730741, 0.26597325539245492, 0.45833053059627044]
    >>> rcorr = LagCorrelationsFromRAP(G0, G1, 3)
    >>> print(rcorr)
    [ 0.01239  0.00274  0.00072]

