butools.map.LagkJointMomentsFromMMAP
====================================

.. currentmodule:: butools.map

.. np:function:: LagkJointMomentsFromMMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Nm = LagkJointMomentsFromMMAP(D, K, L, prec)`
        * - Mathematica:
          - :code:`Nm = LagkJointMomentsFromMMAP[D, K, L, prec]`
        * - Python/Numpy:
          - :code:`Nm = LagkJointMomentsFromMMAP(D, K, L, prec)`

    Returns the lag-L joint moments of a marked Markovian 
    arrival process.

    Parameters
    ----------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the MMAP to check
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
    --------
    For Matlab:

    >>> D0=[-1.78 0.29; 0.07 -0.92];
    >>> D1=[0.15 0.49; 0.23 0.36];
    >>> D2=[0.11 0.2; 0.01 0];
    >>> D3=[0.14 0.4; 0.11 0.14];
    >>> Nm = LagkJointMomentsFromMMAP({D0,D1,D2,D3},4,1);
    >>> Nm{1}
          0.60207      0.60501       1.2755       4.1438
          0.62913      0.62913       1.3226       4.2901
           1.3561       1.3524       2.8387       9.1998
           4.4576       4.4395       9.3105        30.16
    >>> Nm{2}
         0.080053     0.078372      0.16268      0.52401
          0.06033     0.058276      0.11997      0.38467
          0.10244     0.097662       0.1994      0.63637
          0.28923      0.27293       0.5536       1.7601
    >>> Nm{3}
          0.31788      0.31729      0.66629       2.1599
          0.31121      0.30821      0.64424       2.0831
            0.646      0.63672       1.3271       4.2844
           2.0808       2.0455       4.2565        13.73

