butools.ph.AcyclicPHFromME
==========================

.. currentmodule:: butools.ph

.. np:function:: AcyclicPHFromME

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[beta, B] = AcyclicPHFromME(alpha, A, maxSize, precision)`
        * - Mathematica:
          - :code:`{beta, B} = AcyclicPHFromME[alpha, A, maxSize, precision]`
        * - Python/Numpy:
          - :code:`beta, B = AcyclicPHFromME(alpha, A, maxSize, precision)`

    Transforms an arbitrary matrix-exponential representation
    to an acyclic phase-type representation. (see [1]_).
    
    Parameters
    ----------
    alpha : matrix, shape (1,N)
        Initial vector of the distribution
    A : matrix, shape (N,N)
        Matrix parameter of the distribution
    maxSize : int, optional
        The maximum number of phases for the result.
        The default value is 100.
    precision : double, optional
        Vector and matrix entries smaller than the precision
        are considered to be zeros. The default value is 1e-14.
    
    Returns
    -------
    beta : matrix, shape (1,M)
        The initial probability vector of the Markovian 
        acyclic representation
    B : matrix, shape (M,M)
        Transient generator matrix of the Markovian 
        acyclic representation
    
    Notes
    -----
    Raises an error if no Markovian acyclic representation
    has been found.
    
    References
    ----------
    .. [1]  Mocanu, S., Commault, C.: "Sparse representations of
            phase-type distributions," Stoch. Models 15, 759-778 
            (1999)

    Examples
    --------
    For Matlab:

    >>> a=[-0.4 1.4 0];
    >>> A=[-4 1 1; 0 -2 1; 1 0 -8];
    >>> [b,B]=AcyclicPHFromME(a,A);
    >>> b
      0.55273       0.3741     0.073173
    >>> B
      -1.9145       1.9145            0
            0      -3.8858       3.8858
            0            0      -8.1997
    >>> ma=MomentsFromME(a,A,5)
      0.64918      0.73131       1.1825       2.5062       6.5898   
    >>> mb=MomentsFromME(b,B,5);
      0.64918      0.73131       1.1825       2.5062       6.5898
    
    For Python/Numpy:

    >>> a=ml.matrix([[-0.4, 1.4, 0]])
    >>> A=ml.matrix([[-4, 1, 1],[0, -2, 1],[1, 0, -8]])
    >>> b,B=AcyclicPHFromME(a,A)
    >>> print(b)
    [[ 0.55272626  0.37410038  0.07317336]]
    >>> print(B)
    [[-1.91446829  1.91446829  0.        ]
     [ 0.         -3.88582674  3.88582674]
     [ 0.          0.         -8.19970497]]
    >>> print(MomentsFromME(a,A,5))
    [0.64918032786885238, 0.73130878796022558, 1.1825377454500594, 2.5062091064024199, 6.5897616854469261]
    >>> print(MomentsFromME(b,B,5))
    [0.64918032649380331, 0.73130878479925954, 1.182537737569604, 2.5062090835980859, 6.5897616090589066]

