butools.dmap.DMAP2FromMoments
=============================

.. currentmodule:: butools.dmap

.. np:function:: DMAP2FromMoments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[D0, D1] = DMAP2FromMoments(moms, corr1)`
        * - Mathematica:
          - :code:`{D0, D1} = DMAP2FromMoments[moms, corr1]`
        * - Python/Numpy:
          - :code:`D0, D1 = DMAP2FromMoments(moms, corr1)`

    Returns a discrete MAP(2) which has the same 3 marginal
    moments and lag-1 autocorrelation as given.

    Parameters
    ----------
    moms : vector, length(3)
        First three marginal moments of the inter-arrival times
    corr1 : double
        The lag-1 autocorrelation of the inter-arrival times

    Returns
    -------
    D0 : matrix, shape (2,2)
        The D0 matrix of the discrete MAP(2)
    D1 : matrix, shape (2,2)
        The D1 matrix of the discrete MAP(2)

    Notes
    -----
    Raises an exception if the moments are not feasible with
    a DMAP(2). This procedure calls :func:`butools.dmap.DRAPFromMoments`
    followed by :func:`butools.dmap.CanonicalFromDMAP2`.
       
    Examples
    --------
    For Matlab:
    
    >>> moms = [5.1536, 46.587, 626.45];
    >>> corr = -0.00080286;
    >>> [D0,D1]=DMAP2FromMoments(moms,corr);
    >>> D0
              0.3         0.65
          0.61538            0
    >>> D1
             0.05            0
          0.24462         0.14
    >>> MarginalMomentsFromDMAP(D0, D1, 3)
           5.1536       46.587       626.45
    >>> LagCorrelationsFromDMAP(D0, D1, 1)
      -0.00080286

