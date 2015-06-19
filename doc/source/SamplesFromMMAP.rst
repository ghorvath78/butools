butools.map.SamplesFromMMAP
===========================

.. currentmodule:: butools.map

.. np:function:: SamplesFromMMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`x = SamplesFromMMAP(D, K, prec)`
        * - Mathematica:
          - :code:`x = SamplesFromMMAP[D, K, prec]`
        * - Python/Numpy:
          - :code:`x = SamplesFromMMAP(D, K, prec)`

    Generates random samples from a marked Markovian 
    arrival process.

    Parameters
    ----------
    D : list of matrices of shape(M,M), length(N)
        The D0...DN matrices of the MMAP
    K : integer
        The number of samples to generate.
    prec : double, optional
        Numerical precision to check if the input MMAP is
        valid. The default value is 1e-14.

    Returns
    -------
    x : matrix, shape(K,2)
        The random samples. Each row consists of two 
        columns: the inter-arrival time and the type of the
        arrival.        

    Examples
    ========
    For Matlab:

    >>> D0 = [-1.78, 0.29; 0.07, -0.92];
    >>> D1 = [0.15, 0.49; 0.23, 0.36];
    >>> D2 = [0.11, 0.2; 0.01, 0];
    >>> D3 = [0.14, 0.4; 0.11, 0.14];
    >>> D = {D0,D1,D2,D3};
    >>> x = SamplesFromMMAP(D,10000);
    >>> mm = MarginalMomentsFromMMAP(D,3);
    >>> disp(mm);
           1.0007       2.1045       6.8277

    For Mathematica:

    >>> D0 = {{-1.78, 0.29},{0.07, -0.92}};
    >>> D1 = {{0.15, 0.49},{0.23, 0.36}};
    >>> D2 = {{0.11, 0.2},{0.01, 0}};
    >>> D3 = {{0.14, 0.4},{0.11, 0.14}};
    >>> D = {D0,D1,D2,D3};
    >>> x = SamplesFromMMAP[D,10000];
    >>> mm = MarginalMomentsFromMMAP[D,3];
    >>> Print[mm];
    {1.000667111407605, 2.1044966311760755, 6.827688149434602}

    For Python/Numpy:

    >>> D0 = ml.matrix([[-1.78, 0.29],[0.07, -0.92]])
    >>> D1 = ml.matrix([[0.15, 0.49],[0.23, 0.36]])
    >>> D2 = ml.matrix([[0.11, 0.2],[0.01, 0]])
    >>> D3 = ml.matrix([[0.14, 0.4],[0.11, 0.14]])
    >>> D = [D0,D1,D2,D3]
    >>> x = SamplesFromMMAP(D,10000)
    >>> mm = MarginalMomentsFromMMAP(D,3)
    >>> print(mm)
    [1.0006671114076049, 2.1044966311760755, 6.8276881494346]

