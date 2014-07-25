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

