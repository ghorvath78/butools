butools.trace.PdfFromTrace
==========================

.. currentmodule:: butools.trace

.. np:function:: PdfFromTrace

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[x, y] = PdfFromTrace(trace, intBounds)`
        * - Mathematica:
          - :code:`{x, y} = PdfFromTrace[trace, intBounds]`
        * - Python/Numpy:
          - :code:`x, y = PdfFromTrace(trace, intBounds)`

    Returns the empirical density function of a trace.

    Parameters
    ----------
    trace : vector of doubles
        The trace data
    intBounds : vector of doubles
        The array of interval boundaries. The pdf is the
        number of samples falling into an interval divided
        by the interval length.

    Returns
    -------
    x : vector of doubles
        The center of the intervals (the points where the 
        empirical pdf is calculated)
    y : vector of doubles
        The values of the empirical pdf at the given points

    Examples
    ========
    For Matlab:

    >>> D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];
    >>> D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];
    >>> x = (0.0:0.01:0.5);
    >>> tr = SamplesFromMAP(D0, D1, 1000000);
    >>> [x, y] = PdfFromTrace(tr, x);
    >>> [a, A] = MarginalDistributionFromMAP(D0, D1);
    >>> [xm, ym] = IntervalPdfFromPH(a, A, x);
    >>> plot(x, y, xm, ym)

    For Mathematica:

    
    For Python/Numpy:

    >>> D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])
    >>> D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])
    >>> x = np.arange(0.0,0.51,0.01)
    >>> tr = SamplesFromMAP(D0, D1, 1000000)
    >>> x, y = PdfFromTrace(tr, x)
    >>> a, A = MarginalDistributionFromMAP(D0, D1)
    >>> xm, ym = IntervalPdfFromPH(a, A, x)
    >>> plt.plot(x, y, xm, ym)

