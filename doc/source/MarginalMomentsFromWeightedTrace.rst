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
    --------
    For Matlab:
    
    >>> wtr=[0.12; 1.23; 0.546; 0.6765; 1.34; 2.34];
    >>> wei=[12; 1; 34; 23; 8; 2];
    >>> MarginalMomentsFromWeightedTrace(wtr,wei,3)
          0.65242       0.5958      0.74264

    For Python/Numpy:
    
    >>> wtr=np.array([0.12, 1.23, 0.546, 0.6765, 1.34, 2.34])
    >>> wei=np.array([12, 1, 34, 23, 8, 2])
    >>> print(MarginalMomentsFromWeightedTrace(wtr, wei, 3))
    [0.65241875000000005, 0.59579557187499999, 0.74264135759843752]


