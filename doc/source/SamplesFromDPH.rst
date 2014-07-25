butools.dph.SamplesFromDPH
==========================

.. currentmodule:: butools.dph

.. np:function:: SamplesFromDPH

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`x = SamplesFromDPH(alpha, A, K, prec)`
        * - Mathematica:
          - :code:`x = SamplesFromDPH[alpha, A, K, prec]`
        * - Python/Numpy:
          - :code:`x = SamplesFromDPH(alpha, A, K, prec)`

    Generates random samples from a discrete phase-type 
    distribution.

    Parameters
    ----------
    alpha : matrix, shape (1,M)
        The initial probability vector of the discrete phase-
        type distribution.
    A : matrix, shape (M,M)
        The transition probability  matrix of the discrete phase-
        type distribution.
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

    >>> a=[0.76 0 0.24];
    >>> A=[0.34 0.66 0; 0.79 0.05 0.07; 0.26 0.73 0.01];
    >>> x = SamplesFromDPH(a,A,10000000);
    >>> MomentsFromDPH(a,A,5)  
       26.995         1398   1.0853e+05   1.1233e+07   1.4533e+09
    >>> MarginalMomentsFromTrace(x,5)
        26.99         1398    1.086e+05    1.126e+07    1.463e+09
    
