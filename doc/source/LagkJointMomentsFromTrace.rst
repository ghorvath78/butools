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
    --------
    For Matlab:
    
    >>> D0 = [-18 1 4; 2 -18 7; 1 3 -32];
    >>> D1 = [12 1 0; 1 8 0; 2 1 25]; 
    >>> tr = SamplesFromMAP(D0,D1,1000000);
    >>> LagkJointMomentsFromTrace(tr,3,2)
                1      0.05418    0.0064348      0.00121
         0.054181     0.003017   0.00036459   6.9365e-05
        0.0064348   0.00036392   4.4382e-05   8.4913e-06
          0.00121   6.9079e-05   8.4828e-06   1.6295e-06
    >>> LagkJointMomentsFromMAP(D0,D1,3,2)
                1     0.054124     0.006423    0.0012051
         0.054124    0.0030061   0.00036221   6.8492e-05
         0.006423    0.0003622   4.4024e-05    8.361e-06
        0.0012051   6.8486e-05   8.3605e-06   1.5913e-06
        

    For Python/Numpy:
    
    >>> D0 = ml.matrix([[-18, 1, 4],[2, -18, 7],[1, 3, -32]])
    >>> D1 = ml.matrix([[12, 1, 0],[1, 8, 0],[2, 1, 25]])
    >>> tr = SamplesFromMAP(D0,D1,1000000)
    >>> print(LagkJointMomentsFromTrace(tr,3,2))
    [[  1.00000000e+00   5.41123946e-02   6.41988295e-03   1.20353966e-03]
     [  5.41124533e-02   3.01077030e-03   3.63307127e-04   6.87763588e-05]
     [  6.41989123e-03   3.63113615e-04   4.41962106e-05   8.38796443e-06]
     [  1.20354094e-03   6.85773090e-05   8.35435654e-06   1.57742546e-06]]
    >>> print(LagkJointMomentsFromMAP(D0,D1,3,2))
    [[  1.00000000e+00   5.41237113e-02   6.42296483e-03   1.20514648e-03]
     [  5.41237113e-02   3.00610342e-03   3.62213432e-04   6.84916891e-05]
     [  6.42296483e-03   3.62201355e-04   4.40235733e-05   8.36100066e-06]
     [  1.20514648e-03   6.84863446e-05   8.36049904e-06   1.59132033e-06]]
    

