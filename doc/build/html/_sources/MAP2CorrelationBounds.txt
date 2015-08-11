butools.map.MAP2CorrelationBounds
=================================

.. currentmodule:: butools.map

.. np:function:: MAP2CorrelationBounds

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[lb, ub] = MAP2CorrelationBounds(moms)`
        * - Mathematica:
          - :code:`{lb, ub} = MAP2CorrelationBounds[moms]`
        * - Python/Numpy:
          - :code:`lb, ub = MAP2CorrelationBounds(moms)`

    Returns the upper and lower correlation bounds for a MAP(2)
    given the three marginal moments.

    !!!TO CHECK!!!

    Parameters
    ----------
    moms : vector, length(3)
        First three marginal moments of the inter-arrival times

    Returns
    -------
    lb : double
        Lower correlation bound
    ub : double
        Upper correlation bound

    References
    ----------
    .. [1] L Bodrog, A Heindl, G Horvath, M Telek, "A Markovian
           Canonical Form of Second-Order Matrix-Exponential 
           Processes," EUROPEAN JOURNAL OF OPERATIONAL RESEARCH
           190:(2) pp. 459-477. (2008)
           
    Examples
    ========
    For Matlab:

    >>> D0 = [-14., 1.; 1., -25.];
    >>> D1 = [6., 7.; 3., 21.];
    >>> moms = MarginalMomentsFromMAP(D0,D1,3);
    >>> disp(moms);
          0.04918    0.0052609   0.00091819
    >>> [lb,ub] = MAP2CorrelationBounds(moms);
    >>> disp(lb);
        -0.030588
    >>> disp(ub);
         0.074506

    For Mathematica:

    >>> D0 = {{-14., 1.},{1., -25.}};
    >>> D1 = {{6., 7.},{3., 21.}};
    >>> moms = MarginalMomentsFromMAP[D0,D1,3];
    >>> Print[moms];
    {0.04918032786885247, 0.005260932876133214, 0.0009181867601560783}
    >>> {lb,ub} = MAP2CorrelationBounds[moms];
    >>> Print[lb];
    -0.030588145972596268
    >>> Print[ub];
    0.0745055540503923

    For Python/Numpy:

    >>> D0 = ml.matrix([[-14., 1.],[1., -25.]])
    >>> D1 = ml.matrix([[6., 7.],[3., 21.]])
    >>> moms = MarginalMomentsFromMAP(D0,D1,3)
    >>> print(moms)
    [0.049180327868852472, 0.005260932876133214, 0.00091818676015607825]
    >>> lb,ub = MAP2CorrelationBounds(moms)
    >>> print(lb)
    -0.0305881459726
    >>> print(ub)
    0.0745055540504

