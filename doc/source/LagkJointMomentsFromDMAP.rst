butools.dmap.LagkJointMomentsFromDMAP
=====================================

.. currentmodule:: butools.dmap

.. np:function:: LagkJointMomentsFromDMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Nm = LagkJointMomentsFromDMAP(D0, D1, K, L, prec)`
        * - Mathematica:
          - :code:`Nm = LagkJointMomentsFromDMAP[D0, D1, K, L, prec]`
        * - Python/Numpy:
          - :code:`Nm = LagkJointMomentsFromDMAP(D0, D1, K, L, prec)`

    Returns the lag-L joint moments of a discrete Markovian 
    arrival process.

    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the discrete Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the discrete Markovian arrival process
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

    >>> D0 = [0, 0.02, 0, 0; 0, 0.17, 0.2, 0.14; 0.16, 0.17, 0.02, 0.18; 0, 0, 0, 0.12];
    >>> D1 = [0, 0.88, 0.1, 0; 0.18, 0.07, 0.14, 0.1; 0.13, 0.15, 0.15, 0.04; 0.31, 0.18, 0.12, 0.27];
    >>> Nm = LagkJointMomentsFromDMAP(D0, D1, 4, 1);
    >>> disp(Nm);
                1       1.4955       2.9542       7.8852       27.282
           1.4955       2.2037       4.2827       11.293       38.822
           2.9542       4.2875       8.1899       21.315       72.753
           7.8852       11.326       21.379       55.129       187.21
           27.282       38.993        73.17       187.82       636.23
    >>> moms = MarginalMomentsFromDMAP(D0, D1, 4);
    >>> disp(moms);
           1.4955       2.9542       7.8852       27.282
    >>> cjm = zeros(1,3);
    >>> for i=1:1:3
    >>>     Nx = LagkJointMomentsFromDMAP(D0, D1, 1, i);
    >>>     cjm(i) = (Nx(2, 2)-moms(1)^2)/(moms(2)-moms(1)^2);
    >>> end
    >>> disp(cjm);
        -0.045859     0.010753   -0.0047996
    >>> corr = LagCorrelationsFromDMAP(D0, D1, 3);
    >>> disp(corr);
        -0.045859     0.010753   -0.0047996

    For Mathematica:

    
    For Python/Numpy:

    >>> D0 = ml.matrix([[0, 0.02, 0, 0],[0, 0.17, 0.2, 0.14],[0.16, 0.17, 0.02, 0.18],[0, 0, 0, 0.12]])
    >>> D1 = ml.matrix([[0, 0.88, 0.1, 0],[0.18, 0.07, 0.14, 0.1],[0.13, 0.15, 0.15, 0.04],[0.31, 0.18, 0.12, 0.27]])
    >>> Nm = LagkJointMomentsFromDMAP(D0, D1, 4, 1)
    >>> print(Nm)
    [[   1.         1.49554    2.95425    7.88523   27.28233]
     [   1.49554    2.20372    4.28267   11.29332   38.82179]
     [   2.95425    4.28749    8.1899    21.31528   72.75329]
     [   7.88523   11.32649   21.37906   55.12909  187.21291]
     [  27.28233   38.99278   73.17047  187.82218  636.23062]]
    >>> moms = MarginalMomentsFromDMAP(D0, D1, 4)
    >>> print(moms)
    [1.4955358592094412, 2.9542479654368474, 7.885226907678561, 27.282328108669493]
    >>> cjm = np.zeros(3)
    >>> for i in range(1,4,1):
    >>>     Nx = LagkJointMomentsFromDMAP(D0, D1, 1, i)
    >>>     cjm[i-1] = (Nx[1, 1]-moms[0]**2)/(moms[1]-moms[0]**2)
    >>> print(cjm)
    [-0.04586  0.01075 -0.0048 ]
    >>> corr = LagCorrelationsFromDMAP(D0, D1, 3)
    >>> print(corr)
    [-0.04586  0.01075 -0.0048 ]

