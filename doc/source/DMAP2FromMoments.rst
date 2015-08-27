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
    ========
    For Matlab:

    >>> D0 = [0.2, 0.7; 0.6, 0.1];
    >>> D1 = [0.09, 0.01; 0.2, 0.1];
    >>> moms = MarginalMomentsFromDMAP(D0, D1, 3);
    >>> disp(moms);
           5.1536       46.587       626.45
    >>> corr = LagCorrelationsFromDMAP(D0, D1, 1);
    >>> disp(corr);
      -0.00080286
    >>> [D0, D1] = DMAP2FromMoments(moms, corr);
    >>> disp(D0);
              0.3         0.65
          0.61538            0
    >>> disp(D1);
             0.05            0
          0.24462         0.14
    >>> rmoms = MarginalMomentsFromDMAP(D0, D1, 3);
    >>> disp(rmoms);
           5.1536       46.587       626.45
    >>> rcorr = LagCorrelationsFromDMAP(D0, D1, 1);
    >>> disp(rcorr);
      -0.00080286

    For Mathematica:

    >>> D0 = {{0.2, 0.7},{0.6, 0.1}};
    >>> D1 = {{0.09, 0.01},{0.2, 0.1}};
    >>> moms = MarginalMomentsFromDMAP[D0, D1, 3];
    >>> Print[moms];
    {5.15358361774744, 46.58703071672353, 626.4505119453922}
    >>> corr = LagCorrelationsFromDMAP[D0, D1, 1][[1]];
    >>> Print[corr];
    -0.0008028615465149905
    >>> {D0, D1} = DMAP2FromMoments[moms, corr];
    >>> Print[D0];
    {{0.29999999999984595, 0.6500000000001445},
     {0.6153846153846679, 0}}
    >>> Print[D1];
    {{0.050000000000009585, 0},
     {0.24461538461532123, 0.14000000000001087}}
    >>> rmoms = MarginalMomentsFromDMAP[D0, D1, 3];
    >>> Print[rmoms];
    {5.153583617747385, 46.58703071672293, 626.4505119453836}
    >>> rcorr = LagCorrelationsFromDMAP[D0, D1, 1][[1]];
    >>> Print[rcorr];
    -0.0008028615465146371

    For Python/Numpy:

    >>> D0 = ml.matrix([[0.2, 0.7],[0.6, 0.1]])
    >>> D1 = ml.matrix([[0.09, 0.01],[0.2, 0.1]])
    >>> moms = MarginalMomentsFromDMAP(D0, D1, 3)
    >>> print(moms)
    [5.1535836177474383, 46.587030716723511, 626.4505119453919]
    >>> corr = LagCorrelationsFromDMAP(D0, D1, 1)[0]
    >>> print(corr)
    -0.000802861546515
    >>> D0, D1 = DMAP2FromMoments(moms, corr)
    >>> print(D0)
    [[ 0.3      0.65   ]
     [ 0.61538  0.     ]]
    >>> print(D1)
    [[ 0.05     0.     ]
     [ 0.24462  0.14   ]]
    >>> rmoms = MarginalMomentsFromDMAP(D0, D1, 3)
    >>> print(rmoms)
    [5.1535836177474366, 46.587030716723469, 626.45051194539042]
    >>> rcorr = LagCorrelationsFromDMAP(D0, D1, 1)[0]
    >>> print(rcorr)
    -0.000802861546514

