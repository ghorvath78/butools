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
    --------
    For Matlab:

    >>> D0=[-5 0 1 1; 1 -8 1 0; 1 0 -4 1; 1 2 3 -9];
    >>> D1=[0 1 0 2; 2 1 3 0; 0 0 1 1; 1 1 0 1];
    >>> x = SamplesFromMAP(D0,D1,10000000);
    >>> MarginalMomentsFromMAP(D0,D1,5)  
      0.34247      0.25054      0.28271      0.42984      0.81999
    >>> MarginalMomentsFromTrace(x,5)
      0.34246      0.25049      0.28269      0.43013       0.8218
   
    For Python/Numpy:
    
    >>> D0=ml.matrix([[-5, 0, 1, 1],[1, -8, 1, 0],[1, 0, -4, 1],[1, 2, 3, -9]])
    >>> D1=ml.matrix([[0, 1, 0, 2],[2, 1, 3, 0],[0, 0, 1, 1],[1, 1, 0, 1]])
    >>> x = SamplesFromMAP(D0,D1,10000000)
    >>> print(MarginalMomentsFromMAP(D0,D1,5))
    [0.34246575342465752, 0.25053639214391815, 0.28270969431684256, 0.42984404959582057, 0.81998555487921787]    
    >>> print(MarginalMomentsFromTrace(x,5))
    [0.34248531978623137, 0.25064420142125021, 0.28288584219023916, 0.42959312230825203, 0.81564189360938366]

