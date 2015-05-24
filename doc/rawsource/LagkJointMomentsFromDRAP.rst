butools.dmap.LagkJointMomentsFromDRAP
=====================================

.. currentmodule:: butools.dmap

.. np:function:: LagkJointMomentsFromDRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Nm = LagkJointMomentsFromDRAP(H0, H1, K, L, prec)`
        * - Mathematica:
          - :code:`Nm = LagkJointMomentsFromDRAP[H0, H1, K, L, prec]`
        * - Python/Numpy:
          - :code:`Nm = LagkJointMomentsFromDRAP(H0, H1, K, L, prec)`

    Returns the lag-L joint moments of a discrete rational 
    arrival process.

    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the discrete rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the discrete rational arrival process
    K : int, optional
        The dimension of the matrix of joint moments to 
        compute. If K=0, the MxM joint moments will be 
        computed. The default value is 0
    L : int, optional
        The lag at which the joint moments are computed.
        The default value is 1
    prec : double, optional
        Numerical precision to check if the input is valid. 
        The default value is 1e-14

    Returns
    -------
    Nm : matrix, shape(K+1,K+1)
        The matrix containing the lag-L joint moments, 
        starting from moment 0.

    Examples
    --------
    For Matlab:
    
    >>> H0=[0 0 0.13; 0 0.6 0.18; 0.31 0.26 0.02];
    >>> H1=[0 1 -0.13; 0 0.18 0.04; 0.03 0.09 0.29];
    >>> LagkJointMomentsFromDRAP(H0,H1,4,1)
                1        3.207       16.898       130.77       1347.1
            3.207        10.38       54.965       426.11       4391.3
           16.898       54.843       290.86       2256.1        23253
           130.77        424.6       2252.4        17472   1.8009e+05
           1347.1       4373.9        23203   1.7998e+05   1.8552e+06
    >>> MarginalMomentsFromDRAP(H0,H1,4)
            3.207       16.898       130.77       1347.1

    For Python/Numpy:
    
    >>> H0=ml.matrix([[0, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    >>> H1=ml.matrix([[0, 1, -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])
    >>> Nm=LagkJointMomentsFromDRAP(H0,H1,4,1)
    >>> print(Nm)
    [[  1.00000000e+00   3.20702367e+00   1.68976367e+01   1.30770546e+02    1.34707439e+03]
     [  3.20702367e+00   1.03795811e+01   5.49648649e+01   4.26107199e+02    4.39132149e+03]
     [  1.68976367e+01   5.48429776e+01   2.90863052e+02   2.25605103e+03    2.32532721e+04]
     [  1.30770546e+02   4.24600661e+02   2.25239031e+03   1.74717611e+04    1.80086149e+05]
     [  1.34707439e+03   4.37393687e+03   2.32028367e+04   1.79984886e+05    1.85515486e+06]]
    >>> moms=MarginalMomentsFromDRAP(H0,H1,4)
    >>> print(moms)
    [3.2070236684078202, 16.897636691953394, 130.77054574356021, 1347.0743893905096]
    
    
