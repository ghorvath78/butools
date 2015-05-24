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
    ========
    For Matlab:

    >>> a = [0.1, 0.9, 0];
    >>> A = [-6.2, 2, 0; 2, -9, 1; 1, 0, -3];
    >>> x = SamplesFromPH(a,A,1000);
    >>> mp = MomentsFromPH(a,A,3);
    >>> disp(mp);
          0.20939      0.10449     0.089092

    For Mathematica:

    >>> a = {0.1, 0.9, 0};
    >>> A = {{-6.2, 2, 0},{2, -9, 1},{1, 0, -3}};
    >>> x = SamplesFromPH[a,A,1000];
    >>> mp = MomentsFromPH[a,A,3];
    >>> Print[mp];
    {0.20938722294654497, 0.10448912014333092, 0.08909150039190732}

    For Python/Numpy:

    >>> a = ml.matrix([[0.1, 0.9, 0]])
    >>> A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    >>> x = SamplesFromPH(a,A,1000)
    >>> mp = MomentsFromPH(a,A,3)
    >>> print(mp)
    [0.20938722294654497, 0.10448912014333091, 0.089091500391907288]

