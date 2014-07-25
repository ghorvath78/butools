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
    --------
    For Matlab:
    
    >>> moms = [1.5009, 3.1945, 10.026, 43.095, 233.63];
    >>> jmoms = [1, 1.50092, 3.19454; 1.50092, 2.24595, 4.76031; 3.19454, 4.77303, 10.0915];
    >>> [H0,H1]=DRAPFromMoments(moms,jmoms);
    >>> H0
        -0.077384       0.5847     -0.46133
         0.013743     0.097148      0.30219
         0.030442   -0.0019765      0.39024
    >>> H1
          0.23936       2.0767       -1.362
          0.23848      -0.5098      0.85823
          0.19583     -0.10497      0.49044
    >>> MarginalMomentsFromDRAP(H0,H1,5)
           1.5009       3.1945       10.026       43.095       233.63
    >>> LagkJointMomentsFromDRAP(H0,H1,2,1)
                1       1.5009       3.1945
           1.5009       2.2459       4.7603
           3.1945        4.773       10.092

