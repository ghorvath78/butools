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

