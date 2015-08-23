butools.trace.MarginalMomentsFromTrace
======================================

.. currentmodule:: butools.trace

.. np:function:: MarginalMomentsFromTrace

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`moms = MarginalMomentsFromTrace(trace, K)`
        * - Mathematica:
          - :code:`moms = MarginalMomentsFromTrace[trace, K]`
        * - Python/Numpy:
          - :code:`moms = MarginalMomentsFromTrace(trace, K)`

    Returns the marginal moments of a trace.

    Parameters
    ----------
    trace : vector of doubles
        The trace data
    K : int
        The number of moments to compute

    Returns
    -------
    moms : vector of doubles
        The (raw) moments of the trace

    Examples
    ========
    For Matlab:

    >>> D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];
    >>> D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];
    >>> tr = SamplesFromMAP(D0, D1, 1000000);
    >>> moms = MarginalMomentsFromTrace(tr, 3);
    >>> disp(moms);
         0.054215    0.0064528    0.0012149
    >>> mmoms = MarginalMomentsFromMAP(D0, D1, 3);
    >>> disp(mmoms);
         0.054124     0.006423    0.0012051

    For Mathematica:

    
    For Python/Numpy:

    >>> D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])
    >>> D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])
    >>> tr = SamplesFromMAP(D0, D1, 1000000)
    >>> moms = MarginalMomentsFromTrace(tr, 3)
    >>> print(moms)
    [0.054151610439848619, 0.0064330556048736685, 0.001211111157792235]
    >>> mmoms = MarginalMomentsFromMAP(D0, D1, 3)
    >>> print(mmoms)
    [0.054123711340206188, 0.0064229648279705417, 0.0012051464807476149]

