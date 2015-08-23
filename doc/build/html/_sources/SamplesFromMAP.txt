butools.map.SamplesFromMAP
==========================

.. currentmodule:: butools.map

.. np:function:: SamplesFromMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`x = SamplesFromMAP(D0, D1, K, prec)`
        * - Mathematica:
          - :code:`x = SamplesFromMAP[D0, D1, K, prec]`
        * - Python/Numpy:
          - :code:`x = SamplesFromMAP(D0, D1, K, prec)`

    Generates random samples from a Markovian arrival 
    process.

    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process
    K : integer
        The number of samples to generate.
    prec : double, optional
        Numerical precision to check if the input Markovian
        arrival process is valid. The default value is 
        1e-14.

    Returns
    -------
    x : vector, length(K)
        The vector of random samples (inter-arrival times).

    Examples
    ========
    For Matlab:

    >>> D0 = [-0.17, 0, 0, 0.07; 0.01, -0.78, 0.03, 0.08; 0.22, 0.17, -1.1, 0.02; 0.04, 0.12, 0, -0.42];
    >>> D1 = [0, 0.06, 0, 0.04; 0.04, 0.19, 0.21, 0.22; 0.22, 0.13, 0.15, 0.19; 0.05, 0, 0.17, 0.04];
    >>> x = SamplesFromMAP(D0, D1, 10000);
    >>> mt = MarginalMomentsFromTrace(x, 3);
    >>> disp(mt);
            3.458       33.968       568.63
    >>> mm = MarginalMomentsFromMAP(D0, D1, 3);
    >>> disp(mm);
           3.4433        34.03       592.08

    For Mathematica:

    
    For Python/Numpy:

    >>> D0 = ml.matrix([[-0.17, 0, 0, 0.07],[0.01, -0.78, 0.03, 0.08],[0.22, 0.17, -1.1, 0.02],[0.04, 0.12, 0, -0.42]])
    >>> D1 = ml.matrix([[0, 0.06, 0, 0.04],[0.04, 0.19, 0.21, 0.22],[0.22, 0.13, 0.15, 0.19],[0.05, 0, 0.17, 0.04]])
    >>> x = SamplesFromMAP(D0, D1, 10000)
    >>> mt = MarginalMomentsFromTrace(x, 3)
    >>> print(mt)
    [3.4013060257956038, 32.681082499256199, 546.94389771923909]
    >>> mm = MarginalMomentsFromMAP(D0, D1, 3)
    >>> print(mm)
    [3.4433473754205295, 34.030353434960716, 592.07698593172893]

