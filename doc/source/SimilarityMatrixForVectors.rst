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
    ========
    For Matlab:

    >>> vecA = [0.0, 0.3, -1.5, 0.0];
    >>> vecB = [1.0, 0.2, 0.0, 1.0];
    >>> vecA = vecA';
    >>> vecB = vecB';
    >>> B = SimilarityMatrixForVectors (vecA, vecB);
    >>> disp(B);
                0       3.3333            0            0
          0.66667      0.66667            0            0
                0            0            0            0
         -0.83333     -0.83333     -0.83333     -0.83333
    >>> err = norm(B*vecA-vecB);
    >>> disp(err);
         0

    For Mathematica:

    >>> vecA = {0.0, 0.3, -1.5, 0.0};
    >>> vecB = {1.0, 0.2, 0.0, 1.0};
    >>> B = SimilarityMatrixForVectors [vecA, vecB];
    >>> Print[B];
    {{0., 3.3333333333333335, 0., 0.},
     {0.6666666666666667, 0.6666666666666667, 0., 0.},
     {0., 0., 0., 0.},
     {-0.8333333333333334, -0.8333333333333334, -0.8333333333333334, -0.8333333333333334}}
    >>> err = Norm[B.vecA-vecB];
    >>> Print[err];
    0.

    For Python/Numpy:

    >>> vecA = ml.matrix([[0.0, 0.3, -1.5, 0.0]])
    >>> vecB = ml.matrix([[1.0, 0.2, 0.0, 1.0]])
    >>> vecA = vecA.T
    >>> vecB = vecB.T
    >>> B = SimilarityMatrixForVectors (vecA, vecB)
    >>> print(B)
    [[ 0.       3.33333  0.       0.     ]
     [ 0.66667  0.66667  0.       0.     ]
     [ 0.       0.       0.       0.     ]
     [-0.83333 -0.83333 -0.83333 -0.83333]]
    >>> err = la.norm(B*vecA-vecB)
    >>> print(err)
    0.0

