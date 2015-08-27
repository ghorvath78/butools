butools.map.MAP2FromMoments
===========================

.. currentmodule:: butools.map

.. np:function:: MAP2FromMoments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[D0, D1] = MAP2FromMoments(moms, corr1)`
        * - Mathematica:
          - :code:`{D0, D1} = MAP2FromMoments[moms, corr1]`
        * - Python/Numpy:
          - :code:`D0, D1 = MAP2FromMoments(moms, corr1)`

    Returns a MAP(2) which has the same 3 marginal moments 
    and lag-1 autocorrelation as given.

    Parameters
    ----------
    moms : vector, length(3)
        First three marginal moments of the inter-arrival times
    corr1 : double
        The lag-1 autocorrelation of the inter-arrival times

    Returns
    -------
    D0 : matrix, shape (2,2)
        The D0 matrix of the MAP(2)
    D1 : matrix, shape (2,2)
        The D1 matrix of the MAP(2)

    Raises an exception if the moments are not feasible with
    a MAP(2).
    
    Notes
    -----
    The result is always a valid MAP(2) as long as the input
    moments can be realized by a PH(2) (can be checked with 
    :func:`butools.ph.APH2ndMomentLowerBound`, 
    :func:`butools.ph.APH3rdMomentLowerBound`, 
    :func:`butools.ph.APH3rdMomentUpperBound` with n=2) and the 
    correlation falls between the feasible lower and upper 
    bound (check by :func:`MAP2CorrelationBounds`).

    References
    ----------
    .. [1] L Bodrog, A Heindl, G Horvath, M Telek, "A Markovian
           Canonical Form of Second-Order Matrix-Exponential 
           Processes," EUROPEAN JOURNAL OF OPERATIONAL RESEARCH
           190:(2) pp. 459-477. (2008)
        
    Examples
    ========
    For Matlab:

    >>> D0 = [-14., 1.; 1, -25.];
    >>> D1 = [6., 7.; 3., 21.];
    >>> moms = MarginalMomentsFromMAP(D0, D1, 3);
    >>> disp(moms);
          0.04918    0.0052609   0.00091819
    >>> corr = LagCorrelationsFromMAP(D0, D1, 1);
    >>> disp(corr);
         0.022416
    >>> [D0, D1] = MAP2FromMoments(moms, corr);
    >>> disp(D0);
           -13.91        9.199
                0       -25.09
    >>> disp(D1);
           4.7108            0
            2.801       22.289
    >>> rmoms = MarginalMomentsFromMAP(D0, D1, 3);
    >>> disp(rmoms);
          0.04918    0.0052609   0.00091819
    >>> rcorr = LagCorrelationsFromMAP(D0, D1, 1);
    >>> disp(rcorr);
         0.022416

    For Mathematica:

    >>> D0 = {{-14., 1.},{1, -25.}};
    >>> D1 = {{6., 7.},{3., 21.}};
    >>> moms = MarginalMomentsFromMAP[D0, D1, 3];
    >>> Print[moms];
    {0.04918032786885247, 0.005260932876133214, 0.0009181867601560783}
    >>> corr = LagCorrelationsFromMAP[D0, D1, 1][[1]];
    >>> Print[corr];
    0.02241571110398602
    >>> {D0, D1} = MAP2FromMoments[moms, corr];
    >>> Print[D0];
    {{-13.909830056250456, 9.199027971874015},
     {0, -25.090169943749302}}
    >>> Print[D1];
    {{4.710802084376442, 0},
     {2.8009720281259014, 22.2891979156234}}
    >>> rmoms = MarginalMomentsFromMAP[D0, D1, 3];
    >>> Print[rmoms];
    {0.04918032786885251, 0.005260932876133218, 0.0009181867601560789}
    >>> rcorr = LagCorrelationsFromMAP[D0, D1, 1][[1]];
    >>> Print[rcorr];
    0.022415711103985703

    For Python/Numpy:

    >>> D0 = ml.matrix([[-14., 1.],[1, -25.]])
    >>> D1 = ml.matrix([[6., 7.],[3., 21.]])
    >>> moms = MarginalMomentsFromMAP(D0, D1, 3)
    >>> print(moms)
    [0.049180327868852472, 0.005260932876133214, 0.00091818676015607825]
    >>> corr = LagCorrelationsFromMAP(D0, D1, 1)[0]
    >>> print(corr)
    0.022415711104
    >>> D0, D1 = MAP2FromMoments(moms, corr)
    >>> print(D0)
    [[-13.90983   9.19903]
     [  0.      -25.09017]]
    >>> print(D1)
    [[  4.7108    0.     ]
     [  2.80097  22.2892 ]]
    >>> rmoms = MarginalMomentsFromMAP(D0, D1, 3)
    >>> print(rmoms)
    [0.049180327868852479, 0.0052609328761332123, 0.00091818676015607728]
    >>> rcorr = LagCorrelationsFromMAP(D0, D1, 1)[0]
    >>> print(rcorr)
    0.022415711104

