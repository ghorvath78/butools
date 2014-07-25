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
    --------
    For Matlab:
    
    >>> moms = [0.04918, 0.0052609, 0.00091819];
    >>> corr = 0.022416;
    >>> [D0,D1]=MAP2FromMoments(moms,corr);
    >>> D0
           -13.91        9.199
                0       -25.09
    >>> D1
           4.7108            0
            2.801       22.289
    >>> MarginalMomentsFromMAP(D0, D1, 3)
          0.04918    0.0052609   0.00091819
    >>> LagCorrelationsFromMAP(D0, D1, 1)
         0.022416    

