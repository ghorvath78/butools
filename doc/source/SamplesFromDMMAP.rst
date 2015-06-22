butools.dmap.SamplesFromDMMAP
=============================

.. currentmodule:: butools.dmap

.. np:function:: SamplesFromDMMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`x = SamplesFromDMMAP(D, K, prec)`
        * - Mathematica:
          - :code:`x = SamplesFromDMMAP[D, K, prec]`
        * - Python/Numpy:
          - :code:`x = SamplesFromDMMAP(D, K, prec)`

    Generates random samples from a discrete marked 
    Markovian arrival process.

    Parameters
    ----------
    D : list of matrices of shape(M,M), length(N)
        The D0...DN matrices of the DMMAP
    K : integer
        The number of samples to generate.
    prec : double, optional
        Numerical precision to check if the input DMMAP is
        valid. The default value is 1e-14.

    Returns
    -------
    x : matrix, shape(K,2)
        The random samples. Each row consists of two 
        columns: the (discrete) inter-arrival time and the
        type of the arrival.        

    Examples
    ========
    For Matlab:

    >>> D0 = [0.34, 0, 0; 0.06, 0.05, 0.03; 0.11, 0.13, 0];
    >>> D1 = [0.3, 0, 0; 0.16, 0.18, 0.05; 0.15, 0.04, 0.09];
    >>> D2 = [0, 0.01, 0; 0.1, 0.07, 0.08; 0.13, 0.12, 0.13];
    >>> D3 = [0.35, 0, 0; 0, 0.18, 0.04; 0.06, 0.03, 0.01];
    >>> D = {D0,D1,D2,D3};
    >>> x = SamplesFromDMMAP(D,10000);
    >>> mm = MarginalMomentsFromDMMAP(D,3);
    >>> disp(mm);
           1.5037       3.0278       8.4243

    For Mathematica:

    >>> D0 = {{0.34, 0, 0},{0.06, 0.05, 0.03},{0.11, 0.13, 0}};
    >>> D1 = {{0.3, 0, 0},{0.16, 0.18, 0.05},{0.15, 0.04, 0.09}};
    >>> D2 = {{0, 0.01, 0},{0.1, 0.07, 0.08},{0.13, 0.12, 0.13}};
    >>> D3 = {{0.35, 0, 0},{0, 0.18, 0.04},{0.06, 0.03, 0.01}};
    >>> D = {D0,D1,D2,D3};
    >>> x = SamplesFromDMMAP[D,10000];
    >>> mm = MarginalMomentsFromDMMAP[D,3];
    >>> Print[mm];
    {1.503697331491185, 3.027841857350894, 8.424305390832199}

    For Python/Numpy:

    >>> D0 = ml.matrix([[0.34, 0, 0],[0.06, 0.05, 0.03],[0.11, 0.13, 0]])
    >>> D1 = ml.matrix([[0.3, 0, 0],[0.16, 0.18, 0.05],[0.15, 0.04, 0.09]])
    >>> D2 = ml.matrix([[0, 0.01, 0],[0.1, 0.07, 0.08],[0.13, 0.12, 0.13]])
    >>> D3 = ml.matrix([[0.35, 0, 0],[0, 0.18, 0.04],[0.06, 0.03, 0.01]])
    >>> D = [D0,D1,D2,D3]
    >>> x = SamplesFromDMMAP(D,10000)
    >>> mm = MarginalMomentsFromDMMAP(D,3)
    >>> print(mm)
    [1.503697331491185, 3.0278418573508938, 8.424305390832199]

