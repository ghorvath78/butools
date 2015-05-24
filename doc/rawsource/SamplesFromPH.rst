butools.ph.SamplesFromPH
========================

.. currentmodule:: butools.ph

.. np:function:: SamplesFromPH

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`x = SamplesFromPH(alpha, A, K, prec)`
        * - Mathematica:
          - :code:`x = SamplesFromPH[alpha, A, K, prec]`
        * - Python/Numpy:
          - :code:`x = SamplesFromPH(alpha, A, K, prec)`

    Generates random samples from a phase-type distribution.

    Parameters
    ----------
    alpha : matrix, shape (1,M)
        The initial probability vector of the phase-type
        distribution.
    A : matrix, shape (M,M)
        The transient generator matrix of the phase-type
        distribution.
    K : integer
        The number of samples to generate.
    prec : double, optional
        Numerical precision to check if the input phase-type
        distribution is valid. The default value is 1e-14.

    Returns
    -------
    x : vector, length(K)
        The vector of random samples


    Examples
    --------
    For Matlab:

    >>> a = [0.1 0.9 0];
    >>> A = [-6.2 2 0; 2 -9 1; 1 0 -3];
    >>> x = SamplesFromPH(a,A,10000000);
    >>> MomentsFromPH(a,A,5)  
      0.20939      0.10449     0.089092      0.11027      0.17953
    >>> MarginalMomentsFromTrace(x,5)
      0.20942       0.1045     0.089008      0.10977      0.17732
    
    For Python/Numpy:
    
    >>> a=ml.matrix([[0.1, 0.9, 0]])
    >>> A=ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    >>> x = SamplesFromPH(a,A,10000000)
    >>> print(MomentsFromPH(a,A,5))
    [0.20938722294654497, 0.10448912014333091, 0.089091500391907288, 0.11026674096545433, 0.17953027324720897]
    >>> print(MarginalMomentsFromTrace(x,5))
    [0.20947506451314082, 0.10459536790161801, 0.089203808066063053, 0.11031085902290429, 0.17912765854263263]
    
