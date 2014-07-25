butools.dmap.LagkJointMomentsFromDMRAP
======================================

.. currentmodule:: butools.dmap

.. np:function:: LagkJointMomentsFromDMRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Nm = LagkJointMomentsFromDMRAP(H, K, L, prec)`
        * - Mathematica:
          - :code:`Nm = LagkJointMomentsFromDMRAP[H, K, L, prec]`
        * - Python/Numpy:
          - :code:`Nm = LagkJointMomentsFromDMRAP(H, K, L, prec)`

    Returns the lag-L joint moments of a discrete marked 
    rational arrival process.

    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the DMRAP to check
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

    >>> H0=[0.15 0.2 0.18; -0.23 0.17 0.22; 0.19 0.15 0.16];
    >>> H1=[0.01 0.08 0.16; 0.02 0.2 0.07; 0.02 0.15 0.17];
    >>> H2=[0.14 0.07 0.01; 0.19 0.02 0.34; 0.06 0.1 0];
    >>> Nm = LagkJointMomentsFromDMRAP({H0,H1,H2},4,1);
    >>> Nm{1}
          0.48798      0.72584       1.4397       3.9503       14.457
          0.77458       1.1476       2.2651        6.188       22.583
           1.6539       2.4405       4.7925       13.033       47.422
           4.8092       7.0768       13.846       37.534       136.27
           18.218       26.768       52.274       141.45       512.98
    >>> Nm{2}
          0.51202      0.86893       1.9788       6.0092       23.285
          0.82019       1.3944       3.1815       9.6755       37.526
           1.7647        3.005       6.8678       20.914       81.184
           5.1503       8.7784       20.083       61.205       237.71
           19.524       33.289       76.186       232.27       902.29

