butools.trace.PdfFromWeightedTrace
==================================

.. currentmodule:: butools.trace

.. np:function:: PdfFromWeightedTrace

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[x, y] = PdfFromWeightedTrace(trace, weights, intBounds)`
        * - Mathematica:
          - :code:`{x, y} = PdfFromWeightedTrace[trace, weights, intBounds]`
        * - Python/Numpy:
          - :code:`x, y = PdfFromWeightedTrace(trace, weights, intBounds)`

    Returns the empirical density function of a trace 
    consisting of weighted data.

    Parameters
    ----------
    trace : vector of doubles
        The trace data
    weights : vector of doubles
        The weights corresponding to the trace data
    intBounds : vector of doubles
        The array of interval boundaries. The pdf is the 
        number of samples falling into an interval divided by
        the interval length.

    Returns
    -------
    x : vector of doubles
        The center of the intervals (the points where the 
        empirical pdf is calculated)
    y : vector of doubles
        The values of the empirical pdf at the given points

    Examples
    --------
    For Matlab:
    
    >>> wtr=[0.12; 1.23; 0.546; 0.6765; 1.34; 2.34];
    >>> wei=[12; 1; 34; 23; 8; 2];
    >>> [x,y]=PdfFromWeightedTrace(wtr, wei, 0:0.1:3);
    >>> plot(x,y);

    For Python/Numpy:
    
    >>> wtr=np.array([0.12, 1.23, 0.546, 0.6765, 1.34, 2.34])
    >>> wei=np.array([12, 1, 34, 23, 8, 2])
    >>> [x,y]=PdfFromWeightedTrace(wtr, wei, np.arange(0,3.1,0.1))
    >>> plt.plot(x,y)


