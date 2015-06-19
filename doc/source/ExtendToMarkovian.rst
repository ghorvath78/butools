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
    ========
    For Matlab:

    >>> alpha = [0.2, 0.3, 0.5];
    >>> A = [-1, 0, 0; 0, -3, 0.6; 0, -0.6, -3];
    >>> B = TransformToMonocyclic(A);
    >>> disp(B);
               -1            1            0            0
                0      -2.6536       2.6536            0
                0            0      -2.6536       2.6536
                0     0.047227            0      -2.6536
    >>> Cm = SimilarityMatrix(A,B);
    >>> beta = alpha*Cm;
    >>> disp(beta);
         0.045649  -0.00043836    -0.088811       1.0436
    >>> [m,M] = ExtendToMarkovian(beta,B);
    >>> disp(m);
         0.015399     0.010192     0.017621     0.018114    0.0087991      0.10333      0.82655
    >>> disp(M);
               -1            1            0            0            0            0            0
                0      -2.6536       2.6536            0            0            0            0
                0            0      -2.6536       2.6536            0            0            0
                0     0.047227            0      -2.6536       2.6064            0            0
                0            0            0            0      -3.2908       3.2908            0
                0            0            0            0            0      -3.2908       3.2908
                0            0            0            0            0            0      -3.2908
    >>> Cm = SimilarityMatrix(B,M);
    >>> err = norm(B*Cm-Cm*M);
    >>> disp(err);
       5.7254e-15

    For Mathematica:

    >>> alpha = {0.2, 0.3, 0.5};
    >>> A = {{-1, 0, 0},{0, -3, 0.6},{0, -0.6, -3}};
    >>> B = TransformToMonocyclic[A];
    >>> Print[B];
    {{-1., 1., 0, 0},
     {0, -2.6535898384862247, 2.6535898384862247, 0},
     {0, 0, -2.6535898384862247, 2.6535898384862247},
     {0, 0.047227424799192015, 0, -2.6535898384862247}}
    >>> Cm = SimilarityMatrix[A,B];
    >>> beta = alpha.Cm;
    >>> Print[beta];
    {0.0456492140620947, -0.0004383562608946817, -0.08881092881059305, 1.0436000710093931}
    >>> {m,M} = ExtendToMarkovian[beta,B];
    >>> Print[m];
    {0.015398906730148602, 0.010192332530307258, 0.017620759235113224, 0.018114065982724338, 0.008799081502100216, 0.10332742932145772, 0.8265474246981488}
    >>> Print[M];
    {{-1., 1., 0., 0., 0., 0., 0.},
     {0., -2.6535898384862247, 2.6535898384862247, 0., 0., 0., 0.},
     {0., 0., -2.6535898384862247, 2.6535898384862247, 0., 0., 0.},
     {0., 0.047227424799192015, 0., -2.6535898384862247, 2.6063624136870325, 0., 0.},
     {0., 0., 0., 0., -3.290797259447432, 3.290797259447432, 0.},
     {0., 0., 0., 0., 0., -3.290797259447432, 3.290797259447432},
     {0., 0., 0., 0., 0., 0., -3.290797259447432}}
    >>> Cm = SimilarityMatrix[B,M];
    >>> err = Norm[B.Cm-Cm.M];
    >>> Print[err];
    2.7864548696465478*^-15

    For Python/Numpy:

    >>> alpha = ml.matrix([[0.2, 0.3, 0.5]])
    >>> A = ml.matrix([[-1, 0, 0],[0, -3, 0.6],[0, -0.6, -3]])
    >>> B = TransformToMonocyclic(A)
    >>> print(B)
    [[-1.       1.       0.       0.     ]
     [ 0.      -2.65359  2.65359  0.     ]
     [ 0.       0.      -2.65359  2.65359]
     [ 0.       0.04723  0.      -2.65359]]
    >>> Cm = SimilarityMatrix(A,B)
    >>> beta = alpha*Cm
    >>> print(beta)
    [[  4.56492e-02  -4.38356e-04  -8.88109e-02   1.04360e+00]]
    >>> m,M = ExtendToMarkovian(beta,B)
    >>> print(m)
    [[ 0.0154   0.01019  0.01762  0.01811  0.0088   0.10333  0.82655]]
    >>> print(M)
    [[-1.       1.       0.       0.       0.       0.       0.     ]
     [ 0.      -2.65359  2.65359  0.       0.       0.       0.     ]
     [ 0.       0.      -2.65359  2.65359  0.       0.       0.     ]
     [ 0.       0.04723  0.      -2.65359  2.60636  0.       0.     ]
     [ 0.       0.       0.       0.      -3.2908   3.2908   0.     ]
     [ 0.       0.       0.       0.       0.      -3.2908   3.2908 ]
     [ 0.       0.       0.       0.       0.       0.      -3.2908 ]]
    >>> Cm = SimilarityMatrix(B,M)
    >>> err = la.norm(B*Cm-Cm*M)
    >>> print(err)
    4.28823049948e-15

