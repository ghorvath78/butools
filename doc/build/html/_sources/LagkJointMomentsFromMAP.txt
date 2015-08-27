butools.map.LagkJointMomentsFromMAP
===================================

.. currentmodule:: butools.map

.. np:function:: LagkJointMomentsFromMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Nm = LagkJointMomentsFromMAP(D0, D1, K, L, prec)`
        * - Mathematica:
          - :code:`Nm = LagkJointMomentsFromMAP[D0, D1, K, L, prec]`
        * - Python/Numpy:
          - :code:`Nm = LagkJointMomentsFromMAP(D0, D1, K, L, prec)`

    Returns the lag-L joint moments of a Markovian arrival
    process.

    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process
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

    >>> D0 = [-5., 0, 1., 1.; 1., -8., 1., 0; 1., 0, -4., 1.; 1., 2., 3., -9.];
    >>> D1 = [0, 1., 0, 2.; 2., 1., 3., 0; 0, 0, 1., 1.; 1., 1., 0, 1.];
    >>> Nm = LagkJointMomentsFromMAP(D0, D1, 4, 1);
    >>> disp(Nm);
                1      0.34247      0.25054      0.28271      0.42984
          0.34247       0.1173     0.085789     0.096807      0.14721
          0.25054       0.0857     0.062633      0.07066      0.10744
          0.28271     0.096627     0.070589     0.079623      0.12107
          0.42984      0.14686      0.10727      0.12099      0.18396
    >>> moms = MarginalMomentsFromMAP(D0, D1, 4);
    >>> disp(moms);
          0.34247      0.25054      0.28271      0.42984
    >>> cjm = zeros(1,3);
    >>> for i=1:1:3
    >>>     Nx = LagkJointMomentsFromMAP(D0, D1, 1, i);
    >>>     cjm(i) = (Nx(2, 2)-moms(1)^2)/(moms(2)-moms(1)^2);
    >>> end
    >>> disp(cjm);
       0.00012012   0.00086176  -0.00022001
    >>> corr = LagCorrelationsFromMAP(D0, D1, 3);
    >>> disp(corr);
       0.00012012   0.00086176  -0.00022001

    For Mathematica:

    >>> D0 = {{-5., 0, 1., 1.},{1., -8., 1., 0},{1., 0, -4., 1.},{1., 2., 3., -9.}};
    >>> D1 = {{0, 1., 0, 2.},{2., 1., 3., 0},{0, 0, 1., 1.},{1., 1., 0, 1.}};
    >>> Nm = LagkJointMomentsFromMAP[D0, D1, 4, 1];
    >>> Print[Nm];
    {{1., 0.3424657534246575, 0.2505363921439181, 0.2827096943168424, 0.42984404959582045},
     {0.3424657534246575, 0.11729879932812143, 0.08578883767954984, 0.09680718552353199, 0.14720828999251045},
     {0.2505363921439181, 0.08570000543480039, 0.06263282590926178, 0.07065983692223346, 0.10744275082056383},
     {0.28270969431684234, 0.09662651257722407, 0.07058862634724386, 0.07962311566530773, 0.1210669477207951},
     {0.4298440495958204, 0.14686125208953896, 0.10726689466149464, 0.12098747565756454, 0.18395767529024404}}
    >>> moms = MarginalMomentsFromMAP[D0, D1, 4];
    >>> Print[moms];
    {0.3424657534246575, 0.2505363921439181, 0.2827096943168424, 0.42984404959582045}
    >>> cjm = Table[0,{3}];
    >>> Do[
    >>>     Nx = LagkJointMomentsFromMAP[D0, D1, 1, i];
    >>>     cjm[[i]] = (Nx[[2, 2]]-moms[[1]]^2)/(moms[[2]]-moms[[1]]^2);
    >>> , {i,1,3,1}];
    >>> Print[cjm];
    {0.00012012478025432312, 0.0008617649366102103, -0.00022001393374426588}
    >>> corr = LagCorrelationsFromMAP[D0, D1, 3];
    >>> Print[corr];
    {0.00012012478025411484, 0.0008617649366101062, -0.00022001393374437001}

    For Python/Numpy:

    >>> D0 = ml.matrix([[-5., 0, 1., 1.],[1., -8., 1., 0],[1., 0, -4., 1.],[1., 2., 3., -9.]])
    >>> D1 = ml.matrix([[0, 1., 0, 2.],[2., 1., 3., 0],[0, 0, 1., 1.],[1., 1., 0, 1.]])
    >>> Nm = LagkJointMomentsFromMAP(D0, D1, 4, 1)
    >>> print(Nm)
    [[ 1.       0.34247  0.25054  0.28271  0.42984]
     [ 0.34247  0.1173   0.08579  0.09681  0.14721]
     [ 0.25054  0.0857   0.06263  0.07066  0.10744]
     [ 0.28271  0.09663  0.07059  0.07962  0.12107]
     [ 0.42984  0.14686  0.10727  0.12099  0.18396]]
    >>> moms = MarginalMomentsFromMAP(D0, D1, 4)
    >>> print(moms)
    [0.34246575342465752, 0.25053639214391815, 0.28270969431684256, 0.42984404959582057]
    >>> cjm = np.zeros(3)
    >>> for i in range(1,4,1):
    >>>     Nx = LagkJointMomentsFromMAP(D0, D1, 1, i)
    >>>     cjm[i-1] = (Nx[1, 1]-moms[0]**2)/(moms[1]-moms[0]**2)
    >>> print(cjm)
    [ 0.00012  0.00086 -0.00022]
    >>> corr = LagCorrelationsFromMAP(D0, D1, 3)
    >>> print(corr)
    [ 0.00012  0.00086 -0.00022]

