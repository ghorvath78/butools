butools.trace.IATimesFromCummulative
======================================

.. currentmodule:: butools.trace

.. np:function:: IATimesFromCummulative

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`iat = IATimesFromCummulative(trace)`
        * - Mathematica:
          - :code:`iat = IATimesFromCummulative[trace]`
        * - Python/Numpy:
          - :code:`iat = IATimesFromCummulative(trace)`

    Returns the vector of inter-arrival times of a trace
    containing cummulative data.

    Parameters
    ----------
    trace : vector of doubles
        The trace data (cummulative)

    Returns
    -------
    iat : vector of doubles
        The inter-arrival times


