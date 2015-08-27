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
    ========
    For Matlab:

    >>> D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];
    >>> D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];
    >>> tr = SamplesFromMAP(D0, D1, 1000000);
    >>> acf = LagCorrelationsFromTrace(tr, 10);
    >>> disp(acf);
      Columns 1 through 6
         0.043246     0.020421     0.012588    0.0071917    0.0033881    0.0013752
      Columns 7 through 10
        0.0015924  -0.00037795    0.0011809   -0.0002111
    >>> macf = LagCorrelationsFromMAP(D0, D1, 10);
    >>> disp(macf);
      Columns 1 through 6
         0.041288     0.021962     0.012085    0.0068325    0.0039434    0.0023103
      Columns 7 through 10
        0.0013679   0.00081587   0.00048904   0.00029412

    For Mathematica:

    >>> D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};
    >>> D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};
    >>> tr = SamplesFromMAP[D0, D1, 1000000];
    >>> acf = LagCorrelationsFromTrace[tr, 10];
    >>> Print[acf];
    {0.04171700392495588, 0.02236325130968886, 0.011816781009753847, 0.006580180898166265, 0.002870033223830658, 0.0012545187131443494, 0.0015889735659457254, -0.000499385277140403, 0.0006390034516353842, 0.000016288615554384205}
    >>> macf = LagCorrelationsFromMAP[D0, D1, 10];
    >>> Print[macf];
    {0.04128838551112271, 0.021962312595973178, 0.012084710346466325, 0.006832495032493742, 0.003943351206134012, 0.0023102617023690345, 0.0013678941146099855, 0.0008158651499616962, 0.0004890400598515401, 0.0002941204987862201}

    For Python/Numpy:

    >>> D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])
    >>> D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])
    >>> tr = SamplesFromMAP(D0, D1, 1000000)
    >>> acf = LagCorrelationsFromTrace(tr, 10)
    >>> print(acf)
    [0.043090455671182491, 0.023571565460491006, 0.014064127972931892, 0.0085436481559293485, 0.0051488473491149676, 0.0033839569528235254, -0.00053777751233521638, 0.00055429070630664428, 0.00056428151275644746, -0.0018268950662940562]
    >>> macf = LagCorrelationsFromMAP(D0, D1, 10)
    >>> print(macf)
    [ 0.04129  0.02196  0.01208  0.00683  0.00394  0.00231  0.00137  0.00082  0.00049  0.00029]

