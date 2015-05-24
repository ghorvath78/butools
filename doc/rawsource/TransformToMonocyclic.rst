butools.reptrans.TransformToMonocyclic
======================================

.. currentmodule:: butools.reptrans

.. np:function:: TransformToMonocyclic

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`B = TransformToMonocyclic(A, maxSize, precision)`
        * - Mathematica:
          - :code:`B = TransformToMonocyclic[A, maxSize, precision]`
        * - Python/Numpy:
          - :code:`B = TransformToMonocyclic(A, maxSize, precision)`

    Transforms an arbitrary matrix to a Markovian monocyclic 
    matrix (see [1]_).
    
    Parameters
    ----------
    A : matrix, shape (N,N)
        Matrix parameter of the initial representation
    maxSize : int, optional
        The maximal order of the resulting Markovian 
        representation. The default value is 100
    precision : double, optional
        Matrix entries smaller than the precision are 
        considered to be zeros. The default value is 1e-14
        
    Returns
    -------
    B : matrix, shape (M,M)
        Transient generator matrix of the Markovian monocyclic
        representation. Note that M>N if there are complex 
        eigenvalues.

    Notes
    -----    
    Raises an exception if no Markovian monocyclic generator 
    has been found up to the given size.
    
    References
    ----------
    .. [1]  Mocanu, S., Commault, C.: "Sparse representations 
            of phase-type distributions," Stoch. Models 15, 
            759-778 (1999)

    Examples
    --------
    For Matlab:
    
    >>> A = [-1,0,0;0,-3,2;0,-2,-3];
    >>> B=TransformToMonocyclic(A);
       -1        1        0        0        0
        0       -3        3        0        0
        0        0       -3        3        0
        0        0        0       -3        3
        0  0.59259        0        0       -3
    >>> C=SimilarityMatrix(A,B);
    >>> A*C
       1.9399e-15     -0.12308     -0.18462     -0.27692     -0.41538
       6.2774e-17       2.0923      -1.2923      -4.7077       2.9077
       6.6361e-16      0.86154       3.1385      -1.9385      -7.0615
    >>> C*B
       1.9399e-15     -0.12308     -0.18462     -0.27692     -0.41538
       1.1658e-16       2.0923      -1.2923      -4.7077       2.9077
       1.4348e-16      0.86154       3.1385      -1.9385      -7.0615

    For Mathematica:
    
    >>> A = {{-1,0,0},{0,-3,2},{0,-2,-3}};
    >>> B=TransformToMonocyclic[A]
    {{-1, 1, 0, 0, 0}, 
     {0, -3, 3, 0, 0}, 
     {0, 0, -3, 3, 0}, 
     {0, 0, 0, -3, 3}, 
     {0, 16/27, 0, 0, -3}}
    >>> Cm=SimilarityMatrix[A,B];
    >>> A.Cm
    {{0, -(8/65), -(12/65), -(18/65), -(27/65)}, 
    {0, 136/65, -(84/65), -(306/65), 189/65}, 
    {0, 56/65, 204/65, -(126/65), -(459/65)}}
    >>> Cm.B
    {{0, -(8/65), -(12/65), -(18/65), -(27/65)}, 
    {0, 136/65, -(84/65), -(306/65), 189/65}, 
    {0, 56/65, 204/65, -(126/65), -(459/65)}}

    For Python/Numpy:
    
    >>> A = ml.matrix([[-1,0,0],[0,-3,2],[0,-2,-3]])
    >>> B=TransformToMonocyclic(A)
    >>> print(B)
    [[-1.          1.          0.          0.          0.        ]
     [ 0.         -3.          3.          0.          0.        ]
     [ 0.          0.         -3.          3.          0.        ]
     [ 0.          0.          0.         -3.          3.        ]
     [ 0.          0.59259259  0.          0.         -3.        ]]
    >>> C=SimilarityMatrix(A,B)
    >>> print(A*C)
    [[ -3.74874698e-16  -1.23076923e-01  -1.84615385e-01  -2.76923077e-01   -4.15384615e-01]
     [  1.33934087e-15   2.09230769e+00  -1.29230769e+00  -4.70769231e+00    2.90769231e+00]
     [  1.49596465e-15   8.61538462e-01   3.13846154e+00  -1.93846154e+00   -7.06153846e+00]]
    >>> print(C*B)
    [[ -3.74874698e-16  -1.23076923e-01  -1.84615385e-01  -2.76923077e-01   -4.15384615e-01]
     [  5.39227071e-16   2.09230769e+00  -1.29230769e+00  -4.70769231e+00    2.90769231e+00]
     [  1.39170171e-16   8.61538462e-01   3.13846154e+00  -1.93846154e+00   -7.06153846e+00]]

