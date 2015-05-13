butools.reptrans.SimilarityMatrixForVectors
===========================================

.. currentmodule:: butools.reptrans

.. np:function:: SimilarityMatrixForVectors

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`B = SimilarityMatrixForVectors(vecA, vecB)`
        * - Mathematica:
          - :code:`B = SimilarityMatrixForVectors[vecA, vecB]`
        * - Python/Numpy:
          - :code:`B = SimilarityMatrixForVectors(vecA, vecB)`
    
    Returns the similarity transformation matrix that converts 
    a given column vector to an other column vector. It works 
    even with zero entries.
    
    Parameters
    ----------
    vecA : column vector, shape(M,1)
        The original column vector
    vecB : column vector, shape(M,1)
        The target column vector
        
    Returns
    -------
    B : matrix, shape(M,M)
        The matrix by which :math:`B\cdot vecA = vecB` holds

    Examples
    --------
    For Matlab:
    
    >>> vecA = [0.0, 0.3, -1.5, 0.0]';
    >>> vecB = [1.0, 0.2, 0.0, 1.0]';
    >>> B = SimilarityMatrixForVectors (vecA, vecB)
                0       3.3333            0            0
          0.66667      0.66667            0            0
                0            0            0            0
         -0.83333     -0.83333     -0.83333     -0.83333
    >>> B*vecA
                1
              0.2
                0
                1

    For Mathematica:

    >>> vecA = {0.0, 0.3, -1.5, 0.0};
    >>> vecB = {1.0, 0.2, 0.0, 1.0};
    >>> B = SimilarityMatrixForVectors [vecA, vecB]
    {{0., 3.33333, 0., 0.}, 
     {0.666667, 0.666667, 0., 0.}, 
     {0., 0., 0.,  0.}, 
     {-0.833333, -0.833333, -0.833333, -0.833333}}
    >> B.vecA
    {1., 0.2, 0., 1.}        

    For Python/Numpy:
        
    >>> vecA = ml.matrix([[0.0], [0.3], [-1.5], [0.0]])
    >>> vecB = ml.matrix([[1.0], [0.2], [0.0], [1.0]])
    >>> B = SimilarityMatrixForVectors (vecA, vecB)
    >>> print(B)
     [[ 0.          3.33333333  0.          0.        ]
     [ 0.66666667  0.66666667  0.          0.        ]
     [ 0.          0.          0.          0.        ]
     [-0.83333333 -0.83333333 -0.83333333 -0.83333333]]   [[ 0.          3.33333333  0.          0.        ]
    >>> print(B*vecA)
    [[ 1. ]
     [ 0.2]
     [ 0. ]
     [ 1. ]]
    
