butools.trace.CdfFromTrace
==========================

.. currentmodule:: butools.trace

.. np:function:: CdfFromTrace

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[x, y] = CdfFromTrace(trace)`
        * - Mathematica:
          - :code:`{x, y} = CdfFromTrace[trace]`
        * - Python/Numpy:
          - :code:`x, y = CdfFromTrace(trace)`

    Returns the empirical distribution function of the trace.
    
    Parameters
    ----------
    trace : vector of doubles
        The trace data
    
    Returns
    -------
    x : vector of doubles
        The points where the empirical cdf is calculated
    y : vector of doubles
        The values of the empirical cdf at the given points

    Examples
    ========
    For Matlab:

    >>> D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];
    >>> D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];
    >>> tr = SamplesFromMAP(D0, D1, 1000000);
    >>> [x, y] = CdfFromTrace(tr);
    >>> plot(x, y)

    For Mathematica:

    >>> D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};
    >>> D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};
    >>> tr = SamplesFromMAP[D0, D1, 1000000];
    >>> {x, y} = CdfFromTrace[tr];
    >>> ListLinePlot[{Transpose[{x, y}]}]

    For Python/Numpy:

    >>> D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])
    >>> D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])
    >>> tr = SamplesFromMAP(D0, D1, 1000000)
    >>> x, y = CdfFromTrace(tr)
    >>> plt.plot(x, y)

