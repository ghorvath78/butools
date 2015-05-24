butools.trace.LagCorrelationsFromTrace
======================================

.. currentmodule:: butools.trace

.. np:function:: LagCorrelationsFromTrace

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`acf = LagCorrelationsFromTrace(trace, K)`
        * - Mathematica:
          - :code:`acf = LagCorrelationsFromTrace[trace, K]`
        * - Python/Numpy:
          - :code:`acf = LagCorrelationsFromTrace(trace, K)`

    Returns the lag-k autocorrelation of a trace.

    Parameters
    ----------
    trace : vector of doubles
        The trace data
    K : int
        The number of lags to compute

    Returns
    -------
    acf : column vector of doubles
        The lag-k autocorrelation function of the trace up to
        lag K

    Examples
    --------
    For Matlab:
    
    >>> D0 = [-18 1 4; 2 -18 7; 1 3 -32];
    >>> D1 = [12 1 0; 1 8 0; 2 1 25]; 
    >>> tr = SamplesFromMAP(D0,D1,1000000);
    >>> LagCorrelationsFromTrace(tr,5)
         0.040416
         0.023269
         0.011575
        0.0062521
        0.0031942
    >>> LagCorrelationsFromMAP(D0,D1,5)
         0.041288
         0.021962
         0.012085
        0.0068325
        0.0039434

    For Python/Numpy:
    
    >>> D0 = ml.matrix([[-18, 1, 4],[2, -18, 7],[1, 3, -32]])
    >>> D1 = ml.matrix([[12, 1, 0],[1, 8, 0],[2, 1, 25]])
    >>> tr = SamplesFromMAP(D0,D1,1000000)
    >>> print(LagCorrelationsFromTrace(tr,5))
    [0.041536192984685381, 0.023657988505708234, 0.012014712684831945, 0.007522779025904155, 0.0068498273549330788]
    >>> print(LagCorrelationsFromMAP(D0,D1,5))
    [0.041288385511122606, 0.021962312595972942, 0.012084710346466083, 0.0068324950324934984, 0.0039433512061337658]

