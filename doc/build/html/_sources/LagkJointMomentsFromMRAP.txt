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
    ========
    For Matlab:

    >>> x = 0.18;
    >>> H0 = [-5., 0.1+x, 0.9, 1.; 1., -8., 0.9, 0.1; 0.9, 0.1, -4., 1.; 1., 2., 3., -9.];
    >>> H1 = [0.1-x, 0.7, 0.1, 0.1; 0.1, 1., 1.8, 0.1; 0.1, 0.1, 0.1, 0.7; 0.7, 0.1, 0.1, 0.1];
    >>> H2 = [0.1, 0.1, 0.1, 1.7; 1.8, 0.1, 1., 0.1; 0.1, 0.1, 0.7, 0.1; 0.1, 1., 0.1, 0.8];
    >>> Nm = LagkJointMomentsFromMRAP({H0, H1, H2}, 3, 2);
    >>> disp(Nm{1});
          0.41974      0.14337       0.1041      0.11625
          0.14138     0.048248     0.035017       0.0391
          0.10186     0.034737     0.025205      0.02814
          0.11338     0.038655     0.028044     0.031308
    >>> disp(Nm{2});
          0.58026      0.19614      0.14173      0.15799
          0.19813     0.066994     0.048418     0.053974
          0.14397     0.048697     0.035199      0.03924
          0.16086     0.054419     0.039338     0.043855

    For Mathematica:

    >>> x = 0.18;
    >>> H0 = {{-5., 0.1+x, 0.9, 1.},{1., -8., 0.9, 0.1},{0.9, 0.1, -4., 1.},{1., 2., 3., -9.}};
    >>> H1 = {{0.1-x, 0.7, 0.1, 0.1},{0.1, 1., 1.8, 0.1},{0.1, 0.1, 0.1, 0.7},{0.7, 0.1, 0.1, 0.1}};
    >>> H2 = {{0.1, 0.1, 0.1, 1.7},{1.8, 0.1, 1., 0.1},{0.1, 0.1, 0.7, 0.1},{0.1, 1., 0.1, 0.8}};
    >>> Nm = LagkJointMomentsFromMRAP[{H0, H1, H2}, 3, 2];
    >>> Print[Nm[[1]]];
    {{0.41974048988209456, 0.1433714202044304, 0.10409540795191471, 0.11624988638743064},
     {0.14138068443967466, 0.04824753953626901, 0.03501746854001877, 0.039099557275548855},
     {0.10185532272464182, 0.034737397063937, 0.02520539838495395, 0.028140494722472847},
     {0.11337711648763446, 0.038655118635854246, 0.028044484995573313, 0.031308489803818416}}
    >>> Print[Nm[[2]]];
    {{0.5802595101179051, 0.19613605742007034, 0.14173016403045077, 0.15798754127308057},
     {0.1981267931848261, 0.06699407161937415, 0.048417900372854594, 0.05397433134001463},
     {0.14397024925772373, 0.048696896031760306, 0.03519900415612581, 0.039240318805296825},
     {0.16086031117287675, 0.05441861692540258, 0.03933760117012349, 0.0438552148324313}}

    For Python/Numpy:

    >>> x = 0.18
    >>> H0 = ml.matrix([[-5., 0.1+x, 0.9, 1.],[1., -8., 0.9, 0.1],[0.9, 0.1, -4., 1.],[1., 2., 3., -9.]])
    >>> H1 = ml.matrix([[0.1-x, 0.7, 0.1, 0.1],[0.1, 1., 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
    >>> H2 = ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1., 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1., 0.1, 0.8]])
    >>> Nm = LagkJointMomentsFromMRAP([H0, H1, H2], 3, 2)
    >>> print(Nm[0])
    [[ 0.41974  0.14337  0.1041   0.11625]
     [ 0.14138  0.04825  0.03502  0.0391 ]
     [ 0.10186  0.03474  0.02521  0.02814]
     [ 0.11338  0.03866  0.02804  0.03131]]
    >>> print(Nm[1])
    [[ 0.58026  0.19614  0.14173  0.15799]
     [ 0.19813  0.06699  0.04842  0.05397]
     [ 0.14397  0.0487   0.0352   0.03924]
     [ 0.16086  0.05442  0.03934  0.04386]]

