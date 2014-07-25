butools.reptrans.SimilarityMatrix
=================================

.. currentmodule:: butools.reptrans

.. np:function:: SimilarityMatrix

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`B = SimilarityMatrix(A1, A2)`
        * - Mathematica:
          - :code:`B = SimilarityMatrix[A1, A2]`
        * - Python/Numpy:
          - :code:`B = SimilarityMatrix(A1, A2)`
    
    Returns the matrix that transforms A1 to A2.

    Parameters
    ----------
    A1 : matrix, shape (N,N)
        The smaller matrix
    A2 : matrix, shape (M,M)
        The larger matrix (M>=N)
    
    Returns
    -------
    B : matrix, shape (N,M)
        The matrix satisfying :math:`A_1\,B = B\,A_2`
        
    Notes
    -----
    For the existence of a (unique) solution the larger 
    matrix has to inherit the eigenvalues of the smaller one.

    Examples
    --------
    For Matlab:
    
    >>> A1 = [0.2 0.8 0; 1.2 -0.4 0.1; -0.2 0.7 0.5];
    >>> A2 = [0.0088058 0.29904 -0.4672 0.24887; -0.034691 -0.0035465 0.45729 0.02035; -0.20158 1.7687 -1.3378 0.82115; -0.090959 2.2638 -2.2997 1.6325];
    >>> B=SimilarityMatrix(A1,A2)
     0.068433      0.70187     -0.74628      0.97597
     0.030288      0.95884      -1.0099       1.0207
      0.07725      0.71553     -0.80571       1.0129
    >>>  A1*B-B*A2
     -1.9151e-15   8.8818e-16  -7.7716e-16   1.3323e-15
      -9.992e-16   5.5511e-16  -8.8818e-16   1.1102e-15
     -1.5613e-15   1.3323e-15  -6.6613e-16   6.6613e-16

