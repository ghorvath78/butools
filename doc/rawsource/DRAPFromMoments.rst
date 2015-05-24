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

    For Python/Numpy:
    
    >>> moms = [1.4955358592094412, 2.9542479654368474, 7.885226907678561, 27.282328108669493, 116.17171481905851]
    >>> Nm = ml.matrix([[ 1., 1.49553586, 2.95424797],[1.49553586, 2.20371824, 4.2826734],[2.95424797, 4.28748775, 8.18989941]])
    >>> [H0,H1]=DRAPFromMoments(moms,Nm)
    >>> print(H0)
    [[ 0.5644739   0.47187846 -0.69474463]
     [-0.5085687  -0.10550993  0.9592122 ]
     [ 0.18477321  0.26120587 -0.13431385]]
    >>> print(H1)
    [[ 2.39937876  1.12430916 -2.86529566]
     [-1.75345896 -0.59009437  2.99841975]
     [ 0.95074243  0.51878771 -0.78119537]]
    >>> print(MarginalMomentsFromDRAP(H0,H1,5,1e-12))
    [1.495535859209453, 2.9542479654368994, 7.885226907678768, 27.282328108670363, 116.17171481906257]
    >>> print(LagkJointMomentsFromDRAP(H0,H1,2,1,1e-12))
    [[ 1.          1.49553586  2.95424797]
     [ 1.49553586  2.20371824  4.2826734 ]
     [ 2.95424797  4.28748775  8.18989941]]
    
