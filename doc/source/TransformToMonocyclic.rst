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
    ========
    For Matlab:

    >>> A = [-1, 0, 0; 0, -3, 2; 0, -2, -3];
    >>> B = TransformToMonocyclic(A);
    >>> disp(B);
               -1            1            0            0            0
                0           -3            3            0            0
                0            0           -3            3            0
                0            0            0           -3            3
                0      0.59259            0            0           -3
    >>> Cm = SimilarityMatrix(A, B);
    >>> err = norm(A*Cm-Cm*B);
    >>> disp(err);
       1.4306e-14

    For Mathematica:

    >>> A = {{-1, 0, 0},{0, -3, 2},{0, -2, -3}};
    >>> B = TransformToMonocyclic[A];
    >>> Print[B];
    {{-1, 1, 0, 0, 0},
     {0, -3, 3, 0, 0},
     {0, 0, -3, 3, 0},
     {0, 0, 0, -3, 3},
     {0, 16/27, 0, 0, -3}}
    >>> Cm = SimilarityMatrix[A, B];
    >>> err = Norm[A.Cm-Cm.B];
    >>> Print[err];
    0

    For Python/Numpy:

    >>> A = ml.matrix([[-1, 0, 0],[0, -3, 2],[0, -2, -3]])
    >>> B = TransformToMonocyclic(A)
    >>> print(B)
    [[-1.       1.       0.       0.       0.     ]
     [ 0.      -3.       3.       0.       0.     ]
     [ 0.       0.      -3.       3.       0.     ]
     [ 0.       0.       0.      -3.       3.     ]
     [ 0.       0.59259  0.       0.      -3.     ]]
    >>> Cm = SimilarityMatrix(A, B)
    >>> err = la.norm(A*Cm-Cm*B)
    >>> print(err)
    1.29225863835e-14

