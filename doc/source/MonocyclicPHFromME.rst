butools.ph.MonocyclicPHFromME
=============================

.. currentmodule:: butools.ph

.. np:function:: MonocyclicPHFromME

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[beta, B] = MonocyclicPHFromME(alpha, A, maxSize, precision)`
        * - Mathematica:
          - :code:`{beta, B} = MonocyclicPHFromME[alpha, A, maxSize, precision]`
        * - Python/Numpy:
          - :code:`beta, B = MonocyclicPHFromME(alpha, A, maxSize, precision)`

    Transforms an arbitrary matrix-exponential representation
    to a Markovian monocyclic representation.

    Parameters
    ----------
    alpha : matrix, shape (1,N)
        Initial vector of the distribution
    A : matrix, shape (N,N)
        Matrix parameter of the distribution
    maxSize : int, optional
        The maximum number of phases for the result. The default
        value is 100.
    precision : double, optional
        Vector and matrix entries smaller than the precision
        are considered to be zeros. The default value is 1e-14.

    Returns
    -------
    beta : matrix, shape (1,M)
        The initial probability vector of the Markovian 
        monocyclic representation
    B : matrix, shape (M,M)
        Transient generator matrix of the Markovian 
        monocyclic representation

    Notes
    -----
    Raises an error if no Markovian monocyclic representation
    has been found.

    References
    ----------
    .. [1]  Mocanu, S., Commault, C.: "Sparse representations of
           phase-type distributions," Stoch. Models 15, 759-778 
           (1999)

    Examples
    --------
    For Matlab:

    >>> a = [0.2, 0.3, 0.5];
    >>> A = [-1,0,0;0,-3,2;0,-2,-3];
    >>> [b,B]=MonocyclicPHFromME(a,A);
    >>> length(b)
      27
    >>> ma=MomentsFromME(a,A,5)
      0.35385      0.41893       1.1552       4.6998       23.838
    >>> mb=(b,B,5);
      0.35385      0.41893       1.1552       4.6998       23.838

    For Python/Numpy:
    
    >>> a = ml.matrix([[0.2, 0.3, 0.5]])
    >>> A = ml.matrix([[-1,0,0],[0,-3,2],[0,-2,-3]])
    >>> b,B=MonocyclicPHFromME(a,A)
    >>> print(b.shape)
    (1, 27)
    >>> print(MomentsFromME(a,A,5))
    [0.35384615384615381, 0.41893491124260357, 1.1552116522530724, 4.6998354399355771, 23.837756165615836]
    >>> print(MomentsFromME(b,B,5))
    [0.35384615384615531, 0.41893491124260573, 1.155211652253076, 4.699835439935578, 23.83775616561579]

