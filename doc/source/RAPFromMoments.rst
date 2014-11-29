butools.map.RAPFromMoments
==========================

.. currentmodule:: butools.map

.. np:function:: RAPFromMoments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[H0, H1] = RAPFromMoments(moms, Nm)`
        * - Mathematica:
          - :code:`{H0, H1} = RAPFromMoments[moms, Nm]`
        * - Python/Numpy:
          - :code:`H0, H1 = RAPFromMoments(moms, Nm)`

    Creates a rational arrival process that has the same 
    marginal and lag-1 joint moments as given (see [1]_).

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
        The H0 matrix of the rational process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational process
    
    Notes
    -----
    There is no guarantee that the returned matrices define
    a valid stochastic process. The joint densities may be
    negative.

    References
    ----------
    .. [1] G Horvath, M Telek, "A minimal representation of 
           Markov arrival processes and a moments matching 
           method," Performance Evaluation 64:(9-12) pp. 
           1153-1168. (2007)       

    Examples
    --------
    For Matlab:
    
    >>> moms = [0.36585, 0.25535, 0.26507, 0.36691,0.63573];
    >>> jmoms = [1, 0.36585, 0.25535; 0.36585, 0.12866, 0.088334; 0.25535, 0.088802, 0.06067];
    >>> [H0,H1]=RAPFromMoments(moms,jmoms);
    >>> H0
          -12.949        36.78      -24.817
          -1.1102      -2.5113      0.91705
         -0.71205      0.68912      -2.7393
    >>> H1
           9.2672      -99.958       91.678
           1.1693      -2.1771       3.7123
          0.65292       3.9994      -1.8901
    >>> MarginalMomentsFromRAP(H0,H1,5,1e-12)
          0.36585      0.25535      0.26507      0.36691      0.63573
    >>> LagkJointMomentsFromRAP(H0,H1,2,1,1e-12)
                1      0.36585      0.25535
          0.36585      0.12866     0.088334
          0.25535     0.088802      0.06067
    
    For Python/Numpy:
    
    >>> moms = [0.36585365853658536, 0.25535027188603043, 0.26507255497329191, 0.36691170692675046, 0.635727559166956]
    >>> jmoms = ml.matrix([[1., 0.36585366, 0.25535027], [0.36585366, 0.12866312, 0.08833353], [0.25535027, 0.08880202, 0.06066984]])
    >>> [H0,H1]=RAPFromMoments(moms,Nm)
    >>> print(H0)
    [[-12.94938582  36.77982874 -24.81743479]
     [ -1.1101688   -2.51134002   0.91705194]
     [ -0.71205342   0.6891178   -2.73927416]]
    >>> print(H1)
    [[  9.26722501 -99.95808207  91.67784892]
     [  1.16930673  -2.17714855   3.71229869]
     [  0.65291549   3.99937077  -1.89007647]]
    >>> print(MarginalMomentsFromRAP(H0,H1,5,1e-12))
    [0.36585365853658869, 0.25535027188603343, 0.26507255497329529, 0.36691170692675534, 0.635727559166965]
    >>> print(LagkJointMomentsFromRAP(H0,H1,2,1,1e-12))
    [[ 1.          0.36585366  0.25535027]
     [ 0.36585366  0.12866312  0.08833353]
     [ 0.25535027  0.08880202  0.06066984]]

