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
    
