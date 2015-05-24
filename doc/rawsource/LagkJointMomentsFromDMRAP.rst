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

    For Python/Numpy:
    
    >>> H0=ml.matrix([[0.15, 0.2, 0.18],[-0.23, 0.17, 0.22],[0.19, 0.15, 0.16]])
    >>> H1=ml.matrix([[0.01, 0.08, 0.16],[0.02, 0.2, 0.07],[0.02, 0.15, 0.17]])
    >>> H2=ml.matrix([[0.14, 0.07, 0.01],[0.19, 0.02, 0.34],[0.06, 0.1, 0]])
    >>> Nm=LagkJointMomentsFromDMRAP((H0,H1,H2),4,1)
    >>> print(Nm[0])
    [[  4.87980583e-01   7.25835871e-01   1.43972212e+00   3.95026929e+00    1.44570918e+01]
     [  7.74577899e-01   1.14761995e+00   2.26512884e+00   6.18804702e+00    2.25827213e+01]
     [  1.65386099e+00   2.44054938e+00   4.79252795e+00   1.30333369e+01    4.74223768e+01]
     [  4.80922115e+00   7.07677725e+00   1.38464726e+01   3.75336542e+01    1.36274902e+02]
     [  1.82184757e+01   2.67683170e+01   5.22739338e+01   1.41452562e+02    5.12983354e+02]]
    >>> print(Nm[1])
    [[  5.12019417e-01   8.68932390e-01   1.97880994e+00   6.00921631e+00    2.32852777e+01]
     [  8.20190363e-01   1.39442577e+00   3.18148347e+00   9.67551633e+00    3.75263207e+01]
     [  1.76467108e+00   3.00500580e+00   6.86781333e+00   2.09141405e+01    8.11839321e+01]
     [  5.15026444e+00   8.77835281e+00   2.00825841e+01   6.12051218e+01    2.37708056e+02]
     [  1.95238938e+01   3.32889221e+01   7.61863672e+01   2.32268292e+02    9.02287071e+02]]
    
