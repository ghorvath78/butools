butools.trace.MarginalMomentsFromWeightedTrace
==============================================

.. currentmodule:: butools.trace

.. np:function:: MarginalMomentsFromWeightedTrace

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`moms = MarginalMomentsFromWeightedTrace(trace, weights, K)`
        * - Mathematica:
          - :code:`moms = MarginalMomentsFromWeightedTrace[trace, weights, K]`
        * - Python/Numpy:
          - :code:`moms = MarginalMomentsFromWeightedTrace(trace, weights, K)`

    Returns the marginal moments of a trace consisting of 
    weighted data.

    Parameters
    ----------
    trace : vector of doubles
        The trace data
    weights : vector of doubles
        The weights corresponding to the trace data
    K : int
        The number of moments to compute

    Returns
    -------
    moms : vector of doubles
        The (raw) moments of the weighted trace

    Examples
    ========
    For Matlab:

    >>> wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34];
    >>> weights = [12., 1., 34., 23., 8., 2.];
    >>> moms = MarginalMomentsFromWeightedTrace(wtrace, weights, 3);
    >>> disp(moms);
          0.65242       0.5958      0.74264

    For Mathematica:

    >>> wtrace = {0.12, 1.23, 0.546, 0.6765, 1.34, 2.34};
    >>> weights = {12., 1., 34., 23., 8., 2.};
    >>> moms = MarginalMomentsFromWeightedTrace[wtrace, weights, 3];
    >>> Print[moms];
    {0.65241875, 0.595795571875, 0.7426413575984375}

    For Python/Numpy:

    >>> wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34]
    >>> weights = [12., 1., 34., 23., 8., 2.]
    >>> moms = MarginalMomentsFromWeightedTrace(wtrace, weights, 3)
    >>> print(moms)
    [0.65241875000000005, 0.59579557187499999, 0.74264135759843752]

