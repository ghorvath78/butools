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
    ========
    For Matlab:

    >>> H0 = [-2., 0, 0; 0, -3., 1.; 0, -1., -2.];
    >>> H1 = [1.8, 0.2, 0; 0.2, 1.8, 0; 0.2, 1.8, 1.];
    >>> Nm = LagkJointMomentsFromRAP(H0,H1,4,1);
    >>> disp(length(Nm));
         5
    >>> moms = MarginalMomentsFromRAP(H0,H1,4);
    >>> disp(moms);
          0.44444      0.38095      0.48299      0.82216
    >>> cjm=zeros(3,1);
    >>> for i=1:3
    >>>     Nx=LagkJointMomentsFromRAP(H0,H1,1,i);
    >>>     cjm(i) = (Nx(2,2)-moms(1)^2) / (moms(2)-moms(1)^2);
    >>> end
    >>>     Nx=LagkJointMomentsFromRAP(H0,H1,1,i);
    >>>     cjm(i) = (Nx(2,2)-moms(1)^2) / (moms(2)-moms(1)^2);
    >>> end
    >>>     Nx=LagkJointMomentsFromRAP(H0,H1,1,i);
    >>>     cjm(i) = (Nx(2,2)-moms(1)^2) / (moms(2)-moms(1)^2);
    >>> end
    >>> disp(cjm);
       -0.0038462
        0.0045604
        0.0058956
    >>> corr = LagCorrelationsFromRAP(H0,H1,3);
    >>> disp(corr);
       -0.0038462
        0.0045604
        0.0058956

    For Mathematica:

    >>> H0 = {{-2., 0, 0},{0, -3., 1.},{0, -1., -2.}};
    >>> H1 = {{1.8, 0.2, 0},{0.2, 1.8, 0},{0.2, 1.8, 1.}};
    >>> Nm = LagkJointMomentsFromRAP[H0,H1,4,1];
    >>> Print[Length[Nm]];
    5
    >>> moms = MarginalMomentsFromRAP[H0,H1,4];
    >>> Print[moms];
    {0.4444444444444444, 0.380952380952381, 0.48299319727891166, 0.8221574344023325}
    >>> cjm={};
    >>> cjm=Table[(LagkJointMomentsFromRAP[H0,H1,1,i][[2,2]]-moms[[1]]^2) / (moms[[2]]-moms[[1]]^2),{i,1,3}];
    >>> Print[cjm];
    {-0.0038461538461536634, 0.004560439560439573, 0.0058956043956042425}
    >>> corr = LagCorrelationsFromRAP[H0,H1,3];
    >>> Print[corr];
    {-0.0038461538461536634, 0.0045604395604397245, 0.005895604395604545}

    For Python/Numpy:

    >>> H0 = ml.matrix([[-2., 0, 0],[0, -3., 1.],[0, -1., -2.]])
    >>> H1 = ml.matrix([[1.8, 0.2, 0],[0.2, 1.8, 0],[0.2, 1.8, 1.]])
    >>> Nm = LagkJointMomentsFromRAP(H0,H1,4,1)
    >>> print(Length(Nm))
    5
    >>> moms = MarginalMomentsFromRAP(H0,H1,4)
    >>> print(moms)
    [0.44444444444444442, 0.38095238095238093, 0.48299319727891149, 0.82215743440233213]
    >>> cjm=np.empty(3)
    >>> for i in range(3):
    >>>     Nx=LagkJointMomentsFromRAP(H0,H1,1,i+1)
    >>>     cjm[i] = (Nx[1,1]-moms[0]**2) / (moms[1]-moms[0]**2)
    >>>     Nx=LagkJointMomentsFromRAP(H0,H1,1,i+1)
    >>>     cjm[i] = (Nx[1,1]-moms[0]**2) / (moms[1]-moms[0]**2)
    >>>     Nx=LagkJointMomentsFromRAP(H0,H1,1,i+1)
    >>>     cjm[i] = (Nx[1,1]-moms[0]**2) / (moms[1]-moms[0]**2)
    >>> print(cjm)
    [-0.00385  0.00456  0.0059 ]
    >>> corr = LagCorrelationsFromRAP(H0,H1,3)
    >>> print(corr)
    [-0.0038461538461539674, 0.0045604395604395744, 0.005895604395604547]

