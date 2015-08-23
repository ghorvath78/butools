butools.dmap.SamplesFromDMAP
============================

.. currentmodule:: butools.dmap

.. np:function:: SamplesFromDMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`x = SamplesFromDMAP(D0, D1, K, prec)`
        * - Mathematica:
          - :code:`x = SamplesFromDMAP[D0, D1, K, prec]`
        * - Python/Numpy:
          - :code:`x = SamplesFromDMAP(D0, D1, K, prec)`

    Generates random samples from a discrete Markovian 
    arrival process.

    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the discrete MAP.
    D1 : matrix, shape (M,M)
        The D1 matrix of the discrete MAP.
    K : integer
        The number of samples to generate.
    prec : double, optional
        Numerical precision to check if the input DMAP is
        valid. The default value is 1e-14.

    Returns
    -------
    x : vector, length(K)
        The vector of random samples (inter-arrival times).

    Examples
    ========
    For Matlab:

    >>> D0 = [0, 0.02, 0, 0; 0, 0.17, 0.2, 0.14; 0.16, 0.17, 0.02, 0.18; 0, 0, 0, 0.12];
    >>> D1 = [0, 0.88, 0.1, 0; 0.18, 0.07, 0.14, 0.1; 0.13, 0.15, 0.15, 0.04; 0.31, 0.18, 0.12, 0.27];
    >>> x = SamplesFromDMAP(D0, D1, 10000);
    >>> mt = MarginalMomentsFromTrace(x, 3);
    >>> disp(mt);
           1.5171       3.0613       8.3949
    >>> mm = MarginalMomentsFromDMAP(D0, D1, 3);
    >>> disp(mm);
           1.4955       2.9542       7.8852

    For Mathematica:

    
    For Python/Numpy:

    >>> D0 = ml.matrix([[0, 0.02, 0, 0],[0, 0.17, 0.2, 0.14],[0.16, 0.17, 0.02, 0.18],[0, 0, 0, 0.12]])
    >>> D1 = ml.matrix([[0, 0.88, 0.1, 0],[0.18, 0.07, 0.14, 0.1],[0.13, 0.15, 0.15, 0.04],[0.31, 0.18, 0.12, 0.27]])
    >>> x = SamplesFromDMAP(D0, D1, 10000)
    >>> mt = MarginalMomentsFromTrace(x, 3)
    >>> print(mt)
    [1.5088999999999999, 3.0240999999999998, 8.1935000000000002]
    >>> mm = MarginalMomentsFromDMAP(D0, D1, 3)
    >>> print(mm)
    [1.4955358592094412, 2.9542479654368474, 7.885226907678561]

