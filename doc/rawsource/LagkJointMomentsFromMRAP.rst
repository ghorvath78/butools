butools.map.LagkJointMomentsFromMRAP
====================================

.. currentmodule:: butools.map

.. np:function:: LagkJointMomentsFromMRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Nm = LagkJointMomentsFromMRAP(H, K, L, prec)`
        * - Mathematica:
          - :code:`Nm = LagkJointMomentsFromMRAP[H, K, L, prec]`
        * - Python/Numpy:
          - :code:`Nm = LagkJointMomentsFromMRAP(H, K, L, prec)`

    Returns the lag-L joint moments of a marked rational 
    arrival process.

    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the MRAP to check
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
    Nm : list/cell of matrices of shape(K+1,K+1), length(N)
        The matrices containing the lag-L joint moments,
        starting from moment 0.

    Examples
    --------
    For Matlab:

    >>> H0=[-5 0.28 0.9 1; 1 -8 0.9 0.1; 0.9 0.1 -4 1; 1 2 3 -9]
    >>> H1=[-0.08 0.7 0.1 0.1; 0.1 1 1.8 0.1; 0.1 0.1 0.1 0.7; 0.7 0.1 0.1 0.1]
    >>> H2=[0.1 0.1 0.1 1.7; 1.8 0.1 1 0.1; 0.1 0.1 0.7 0.1; 0.1 1 0.1 0.8]
    >>> Nm = LagkJointMomentsFromRMAP({H0,H1,H2},4,1);
    >>> Nm{1}
          0.41974      0.14337       0.1041      0.11625
          0.14138     0.048248     0.035017       0.0391
          0.10186     0.034737     0.025205      0.02814
          0.11338     0.038655     0.028044     0.031308
    >>> Nm{2}
          0.58026      0.19614      0.14173      0.15799
          0.19813     0.066994     0.048418     0.053974
          0.14397     0.048697     0.035199      0.03924
          0.16086     0.054419     0.039338     0.043855

    For Python/Numpy:
    
    >>> H0=ml.matrix([[-5, 0.28, 0.9, 1],[1, -8, 0.9, 0.1],[0.9, 0.1, -4, 1],[1, 2, 3, -9]])
    >>> H1=ml.matrix([[-0.08, 0.7, 0.1, 0.1],[0.1, 1, 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
    >>> H2=ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1, 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1, 0.1, 0.8]])
    >>> Nm=LagkJointMomentsFromMRAP((H0,H1,H2),3,1)
    >>> print(Nm[0])
    [[ 0.41974049  0.13785882  0.09865873  0.10965738]
     [ 0.14138068  0.04586476  0.03262009  0.03616431]
     [ 0.10185532  0.03287395  0.02331696  0.02581984]
     [ 0.11337712  0.03652543  0.02588026  0.02864508]]
    >>> print(Nm[1])
    [[ 0.58025951  0.20164866  0.14716684  0.16458005]
     [ 0.19812679  0.06934517  0.05076715  0.05685583]
     [ 0.14397025  0.05046534  0.03697262  0.04142364]
     [ 0.16086031  0.05639756  0.04132497  0.046305  ]]

    
