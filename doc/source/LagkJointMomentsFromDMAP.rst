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
    --------
    For Matlab:
    
    >>> D0=[0 0.02 0 0; 0 0.17 0.2 0.14; 0.16 0.17 0.02 0.18; 0 0 0 0.12];
    >>> D1=[0 0.88 0.1 0; 0.18 0.07 0.14 0.1; 0.13 0.15 0.15 0.04; 0.31 0.18 0.12 0.27];
    >>> LagkJointMomentsFromDMAP(D0,D1,4,1)
                1       1.4955       2.9542       7.8852       27.282
           1.4955       2.2037       4.2827       11.293       38.822
           2.9542       4.2875       8.1899       21.315       72.753
           7.8852       11.326       21.379       55.129       187.21
           27.282       38.993        73.17       187.82       636.23
    >>> MarginalMomentsFromDMAP(D0,D1,4)
           1.4955       2.9542       7.8852       27.282

