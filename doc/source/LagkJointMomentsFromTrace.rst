butools.trace.LagCorrelationsFromTrace
======================================

.. currentmodule:: butools.trace

.. np:function:: LagkJointMomentsFromTrace

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Nm = LagkJointMomentsFromTrace(trace, K, L)`
        * - Mathematica:
          - :code:`Nm = LagkJointMomentsFromTrace[trace, K, L]`
        * - Python/Numpy:
          - :code:`Nm = LagkJointMomentsFromTrace(trace, K, L)`

    Returns the lag-L joint moments of a trace.

    It is computed as :math:`Nm_{i,j}=\frac{1}{N-L}\sum_{k=0}^{N-L} x_k^i x_{k+L}^j`.

    Parameters
    ----------
    trace : vector of doubles
        The trace data
    K : int
        The joint moments are computed up to order K
    L : int, optional
        The lag-L joint moments are computed.
        The default value is 1.

    Returns
    -------
    Nm : matrix, shape (K,K)
        The matrix of joint moments, starting from moment 0

    Examples
    ========
    For Matlab:

    >>> D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];
    >>> D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];
    >>> tr = SamplesFromMAP(D0, D1, 1000000);
    >>> Nm1 = LagkJointMomentsFromTrace(tr, 3, 1);
    >>> disp(Nm1);
                1     0.054231    0.0064413    0.0012077
         0.054231    0.0030911    0.0003777   7.1599e-05
        0.0064413   0.00037728   4.6603e-05   8.8002e-06
        0.0012077   7.1519e-05   8.8317e-06   1.6498e-06
    >>> mNm1 = LagkJointMomentsFromMAP(D0, D1, 3, 1);
    >>> disp(mNm1);
                1     0.054124     0.006423    0.0012051
         0.054124    0.0030736   0.00037516   7.1415e-05
         0.006423   0.00037515   4.6507e-05   8.9222e-06
        0.0012051    7.141e-05   8.9217e-06   1.7182e-06
    >>> Nm2 = LagkJointMomentsFromTrace(tr, 3, 2);
    >>> disp(Nm2);
                1     0.054231    0.0064413    0.0012077
         0.054231    0.0030213   0.00036415   6.8624e-05
        0.0064414   0.00036417   4.4121e-05   8.2941e-06
        0.0012077   6.8666e-05   8.2925e-06   1.5427e-06
    >>> mNm2 = LagkJointMomentsFromMAP(D0, D1, 3, 2);
    >>> disp(mNm2);
                1     0.054124     0.006423    0.0012051
         0.054124    0.0030061   0.00036221   6.8492e-05
         0.006423    0.0003622   4.4024e-05    8.361e-06
        0.0012051   6.8486e-05   8.3605e-06   1.5913e-06

    For Mathematica:

    
    For Python/Numpy:

    >>> D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])
    >>> D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])
    >>> tr = SamplesFromMAP(D0, D1, 1000000)
    >>> Nm1 = LagkJointMomentsFromTrace(tr, 3, 1)
    >>> print(Nm1)
    [[  1.00000e+00   5.41472e-02   6.42593e-03   1.20511e-03]
     [  5.41473e-02   3.07546e-03   3.75237e-04   7.12258e-05]
     [  6.42593e-03   3.75451e-04   4.64806e-05   8.87509e-06]
     [  1.20511e-03   7.14549e-05   8.89829e-06   1.70148e-06]]
    >>> mNm1 = LagkJointMomentsFromMAP(D0, D1, 3, 1)
    >>> print(mNm1)
    [[  1.00000e+00   5.41237e-02   6.42296e-03   1.20515e-03]
     [  5.41237e-02   3.07362e-03   3.75157e-04   7.14155e-05]
     [  6.42296e-03   3.75145e-04   4.65066e-05   8.92218e-06]
     [  1.20515e-03   7.14101e-05   8.92168e-06   1.71822e-06]]
    >>> Nm2 = LagkJointMomentsFromTrace(tr, 3, 2)
    >>> print(Nm2)
    [[  1.00000e+00   5.41472e-02   6.42593e-03   1.20511e-03]
     [  5.41472e-02   3.00173e-03   3.61013e-04   6.79513e-05]
     [  6.42591e-03   3.60964e-04   4.37219e-05   8.21233e-06]
     [  1.20511e-03   6.82216e-05   8.30943e-06   1.55532e-06]]
    >>> mNm2 = LagkJointMomentsFromMAP(D0, D1, 3, 2)
    >>> print(mNm2)
    [[  1.00000e+00   5.41237e-02   6.42296e-03   1.20515e-03]
     [  5.41237e-02   3.00610e-03   3.62213e-04   6.84917e-05]
     [  6.42296e-03   3.62201e-04   4.40236e-05   8.36100e-06]
     [  1.20515e-03   6.84863e-05   8.36050e-06   1.59132e-06]]

