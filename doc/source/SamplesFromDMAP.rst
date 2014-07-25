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
    --------
    For Matlab:

    >>> D0=[0 0.02 0 0; 0 0.17 0.2 0.14; 0.16 0.17 0.02 0.18; 0 0 0 0.12];
    >>> D1=[0 0.88 0.1 0; 0.18 0.07 0.14 0.1; 0.13 0.15 0.15 0.04; 0.31 0.18 0.12 0.27];
    >>> x = SamplesFromDMAP(D0,D1,10000000);
    >>> MarginalMomentsFromDMAP(D0,D1,5)  
       1.4955       2.9542       7.8852       27.282       116.17
    >>> MarginalMomentsFromTrace(x,5)
       1.4958       2.9551        7.888       27.292       116.27
    
