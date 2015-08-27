butools.trace.CdfFromWeightedTrace
==================================

.. currentmodule:: butools.trace

.. np:function:: CdfFromWeightedTrace

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[x, y] = CdfFromWeightedTrace(trace, weights)`
        * - Mathematica:
          - :code:`{x, y} = CdfFromWeightedTrace[trace, weights]`
        * - Python/Numpy:
          - :code:`x, y = CdfFromWeightedTrace(trace, weights)`

    Returns the empirical distribution function of a trace
    consisting of weighted data.
    
    Parameters
    ----------
    trace : vector of doubles
        The trace data
    weights : vector of doubles
        The weights corresponding to the trace data
    
    Returns
    -------
    x : vector of doubles
        The points where the empirical cdf is calculated
    y : vector of doubles
        The values of the empirical cdf at the given points

    Examples
    ========
    For Matlab:

    >>> wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34];
    >>> weights = [12., 1., 34., 23., 8., 2.];
    >>> [x, y] = CdfFromWeightedTrace(wtrace, weights);
    >>> plot(x, y)

    For Mathematica:

    >>> wtrace = {0.12, 1.23, 0.546, 0.6765, 1.34, 2.34};
    >>> weights = {12., 1., 34., 23., 8., 2.};
    >>> {x, y} = CdfFromWeightedTrace[wtrace, weights];
    >>> ListLinePlot[{Transpose[{x, y}]}]

    For Python/Numpy:

    >>> wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34]
    >>> weights = [12., 1., 34., 23., 8., 2.]
    >>> x, y = CdfFromWeightedTrace(wtrace, weights)
    >>> plt.plot(x, y)

