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
    --------
    For Matlab:

    >>> D0=[-1.78 0.29; 0.07 -0.92]
    >>> D1=[0.15 0.49; 0.23 0.36]
    >>> D2=[0.11 0.2; 0.01 0]
    >>> D3=[0.14 0.4; 0.11 0.14]
    >>> D = {D0,D1,D2,D3};
    >>> x = SamplesFromMMAP(D,10000000);
    >>> MarginalMomentsFromMMAP(D,5)  
       1.0007       2.1045       6.8277       29.979       165.81
    >>> MarginalMomentsFromTrace(x(:,1),5)
       1.0014       2.1077       6.8466       30.115       166.94
   
    For Python/Numpy:
    
    >>> D0=ml.matrix([[-1.78, 0.29],[0.07, -0.92]])
    >>> D1=ml.matrix([[0.15, 0.49],[0.23, 0.36]])
    >>> D2=ml.matrix([[0.11, 0.2],[0.01, 0]])
    >>> D3=ml.matrix([[0.14, 0.4],[0.11, 0.14]])
    >>> D = (D0,D1,D2,D3)
    >>> x = SamplesFromMMAP(D,10000000)
    >>> print(MarginalMomentsFromMMAP(D,5))
    [1.0006671114076049, 2.1044966311760755, 6.8276881494346, 29.979022098423791, 165.80908634145123]
    >>> print(MarginalMomentsFromTrace(x[:,0],5))
    [1.0007512647897872, 2.1050907410974111, 6.8318440377472891, 30.016797811190671, 166.22470555982869]

    
