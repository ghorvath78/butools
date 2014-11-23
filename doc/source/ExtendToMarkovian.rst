butools.reptrans.ExtendToMarkovian
===================================

.. currentmodule:: butools.reptrans

.. np:function:: ExtendToMarkovian

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[beta, B] = ExtendToMarkovian(alpha, A, maxSize, precision)`
        * - Mathematica:
          - :code:`{beta, B} = ExtendToMarkovian[alpha, A, maxSize, precision]`
        * - Python/Numpy:
          - :code:`beta, B = ExtendToMarkovian(alpha, A, maxSize, precision)`

    Assume we have an existing monocyclic (or acyclic) 
    representation of a matrix-exponential distribution 
    described by matrix A and vector alpha such that A is 
    Markovian but alpha is not. This procedure appends an 
    appropriate Erlang tail to the representation that makes 
    the result Markovian (both the generator matrix and the 
    initial vector parameter), while keeping the distribution 
    the same. In [1]_ it is proven that this is always 
    possible if the initial (alpha,A) defines a distribuion 
    (non-negative density).
    
    Parameters
    ----------
    alpha : vector, shape (1,M)
        The (non-Markovian) initial vector
    A : matrix, shape (M,M)            
        The (Markovian) transient generator.
    maxSize : int, optional
        The procedure stops if more than maxSize new phases 
        are required. The default value is 100
    precision : double, optional
        The initial vector is considered to be valid if the 
        smallest entry is greater than -precision. The
        default value is 1e-14
    
    Returns
    -------
    beta : vector, shape (1,N)
        The Markovian initial vector (N>=M)
    B : matrix, shape (N,N)
        The Markovian transient generator (N>=M).
    
    References
    ----------
    .. [1]  Mocanu, S., Commault, C.: "Sparse representations 
            of phase-type distributions," Stoch. Models 15, 
            759-778 (1999)

    Examples
    --------
    
    This example transforms a matrix-exponential distribution 
    (gamma,G) to monocyclic representation (alpha,A). The 
    resulting initial vector is, however, non Markovian 
    (observe the negative entries). The ExtendToMarkovian
    function is able to add an Erlang tail such that it gets
    Markovian, too, while keeping the distribution the same.
    
    For Matlab:
    
    >>> gamma = [0.2, 0.3, 0.5];
    >>> G = [-1,0,0;0,-3,0.6;0,-0.6,-3];
    >>> A=TransformToMonocyclic(G);
       -1            1            0            0
        0      -2.6536       2.6536            0
        0            0      -2.6536       2.6536
        0     0.047227            0      -2.6536
    >>> C=SimilarityMatrix(G,A);
    >>> alpha = gamma*C
        0.045649  -0.00043836    -0.088811       1.0436
    >>> [beta,B]=ExtendToMarkovian(alpha,A);
    >>> beta
        0.015399     0.010192     0.017621     0.018114    0.0087991      0.10333      0.82655
    >>> B
       -1            1            0            0            0            0            0
        0      -2.6536       2.6536            0            0            0            0
        0            0      -2.6536       2.6536            0            0            0
        0     0.047227            0      -2.6536       2.6064            0            0
        0            0            0            0      -3.2908       3.2908            0
        0            0            0            0            0      -3.2908       3.2908
        0            0            0            0            0            0      -3.2908
    
    For Python/Numpy:
    
    >>> gamma = ml.matrix([[0.2, 0.3, 0.5]])    
    >>> G = ml.matrix([[-1,0,0],[0,-3,0.6],[0,-0.6,-3]])
    >>> A=TransformToMonocyclic(G)
    >>> print(A)
    [[-1.          1.          0.          0.        ]
     [ 0.         -2.65358984  2.65358984  0.        ]
     [ 0.          0.         -2.65358984  2.65358984]
     [ 0.          0.04722742  0.         -2.65358984]]
    >>> C=SimilarityMatrix(G,A)
    >>> alpha = gamma*C
    >>> print(alpha)
    [[  4.56492145e-02  -4.38356417e-04  -8.88109272e-02   1.04360007e+00]]
    >>> beta, B = ExtendToMarkovian(alpha,A)
    >>> print(beta)
    [[ 0.01539891  0.01019233  0.01762076  0.01811407  0.00879908  0.10332743   0.82654742]]   
    >>> print(B)
    [[-1.          1.          0.          0.          0.          0.          0.        ]
     [ 0.         -2.65358984  2.65358984  0.          0.          0.          0.        ]
     [ 0.          0.         -2.65358984  2.65358984  0.          0.          0.        ]
     [ 0.          0.04722742  0.         -2.65358984  2.60636242  0.          0.        ]
     [ 0.          0.          0.          0.         -3.29079727  3.29079727   0.        ]
     [ 0.          0.          0.          0.          0.         -3.29079727   3.29079727]
     [ 0.          0.          0.          0.          0.          0.  -3.29079727]]
    
