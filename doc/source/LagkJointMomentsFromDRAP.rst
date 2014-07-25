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

