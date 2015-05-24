butools.map.LagkJointMomentsFromRAP
===================================

.. currentmodule:: butools.map

.. np:function:: LagkJointMomentsFromRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Nm = LagkJointMomentsFromRAP(H0, H1, K, L, prec)`
        * - Mathematica:
          - :code:`Nm = LagkJointMomentsFromRAP[H0, H1, K, L, prec]`
        * - Python/Numpy:
          - :code:`Nm = LagkJointMomentsFromRAP(H0, H1, K, L, prec)`

    Returns the lag-L joint moments of a rational arrival
    process.

    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational arrival process
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
    
    >>> H0=[-2 0 0; 0 -3 1; 0 -1 -2];
    >>> H1=[1.8 0.2 0; 0.2 1.8 0; 0.2 1.8 1];
    >>> LagkJointMomentsFromRAP(H0,H1,4,1)
                1      0.44444      0.38095      0.48299      0.82216
          0.44444      0.19683      0.16961      0.21788      0.37801
          0.38095      0.16984      0.14879      0.19574      0.34921
          0.48299      0.21871      0.19643      0.26642      0.49043
          0.82216      0.38066      0.35212       0.4931      0.93488
    >>> MarginalMomentsFromRAP(H0,H1,4)
          0.44444      0.38095      0.48299      0.82216

    For Python/Numpy:
    
    >>> H0=ml.matrix([[-2, 0, 0],[0, -3, 1],[0, -1, -2]])
    >>> H1=ml.matrix([[1.8, 0.2, 0],[0.2, 1.8, 0],[0.2, 1.8, 1]])
    >>> print(LagkJointMomentsFromRAP(H0,H1,4,1))
    [[ 1.          0.44444444  0.38095238  0.4829932   0.82215743]
     [ 0.44444444  0.1968254   0.16961451  0.21788144  0.37800916]
     [ 0.38095238  0.16984127  0.14878523  0.19574483  0.34920569]
     [ 0.4829932   0.21870748  0.1964251   0.26642322  0.49043213]
     [ 0.82215743  0.38066084  0.35212114  0.49309769  0.93488051]]
    >>> print(MarginalMomentsFromRAP(H0,H1,4))
    [0.44444444444444442, 0.38095238095238093, 0.48299319727891149, 0.82215743440233213]

