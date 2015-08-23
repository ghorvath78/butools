butools.dmap.LagkJointMomentsFromDMMAP
======================================

.. currentmodule:: butools.dmap

.. np:function:: LagkJointMomentsFromDMMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Nm = LagkJointMomentsFromDMMAP(D, K, L, prec)`
        * - Mathematica:
          - :code:`Nm = LagkJointMomentsFromDMMAP[D, K, L, prec]`
        * - Python/Numpy:
          - :code:`Nm = LagkJointMomentsFromDMMAP(D, K, L, prec)`

    Returns the lag-L joint moments of a discrete marked 
    Markovian arrival process.

    Parameters
    ----------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the DMMAP to check
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
    Nm : list/cell of matrices of shape(K+1,K+1), length(L)
        The matrices containing the lag-L joint moments,
        starting from moment 0.

    Examples
    ========
    For Matlab:

    >>> D0 = [0.34, 0, 0; 0.06, 0.05, 0.03; 0.11, 0.13, 0];
    >>> D1 = [0.3, 0, 0; 0.16, 0.18, 0.05; 0.15, 0.04, 0.09];
    >>> D2 = [0, 0.01, 0; 0.1, 0.07, 0.08; 0.13, 0.12, 0.13];
    >>> D3 = [0.35, 0, 0; 0, 0.18, 0.04; 0.06, 0.03, 0.01];
    >>> Nm = LagkJointMomentsFromDMMAP({D0, D1, D2, D3}, 3, 1);
    >>> disp(Nm{1});
          0.45395      0.68525       1.3856       3.8671
          0.68283       1.0318       2.0887       5.8339
           1.3755       2.0807       4.2171       11.789
            3.828       5.7954       11.756       32.887
    >>> disp(Nm{2});
         0.026281      0.03323     0.053055      0.11925
         0.035051     0.043866     0.068917      0.15222
         0.060653      0.07477      0.11464      0.24631
           0.1482      0.17996      0.26901      0.56074
    >>> disp(Nm{3});
          0.51977      0.78522       1.5891        4.438
          0.78582       1.1881       2.4067       6.7254
           1.5917       2.4087       4.8838       13.657
           4.4481       6.7354       13.666       38.235

    For Mathematica:

    
    For Python/Numpy:

    >>> D0 = ml.matrix([[0.34, 0, 0],[0.06, 0.05, 0.03],[0.11, 0.13, 0]])
    >>> D1 = ml.matrix([[0.3, 0, 0],[0.16, 0.18, 0.05],[0.15, 0.04, 0.09]])
    >>> D2 = ml.matrix([[0, 0.01, 0],[0.1, 0.07, 0.08],[0.13, 0.12, 0.13]])
    >>> D3 = ml.matrix([[0.35, 0, 0],[0, 0.18, 0.04],[0.06, 0.03, 0.01]])
    >>> Nm = LagkJointMomentsFromDMMAP([D0, D1, D2, D3], 3, 1)
    >>> print(Nm[0])
    [[  0.45395   0.68525   1.38564   3.86705]
     [  0.68283   1.03179   2.08868   5.83387]
     [  1.37545   2.08069   4.21706  11.78914]
     [  3.82802   5.79545  11.75638  32.88737]]
    >>> print(Nm[1])
    [[ 0.02628  0.03323  0.05306  0.11925]
     [ 0.03505  0.04387  0.06892  0.15222]
     [ 0.06065  0.07477  0.11464  0.24631]
     [ 0.1482   0.17996  0.26901  0.56074]]
    >>> print(Nm[2])
    [[  0.51977   0.78522   1.58915   4.438  ]
     [  0.78582   1.1881    2.40665   6.72537]
     [  1.59173   2.4087    4.88378  13.65717]
     [  4.44809   6.73539  13.66584  38.23481]]

