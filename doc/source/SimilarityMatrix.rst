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
    ========
    For Matlab:

    >>> A1m = [0.2, 0.8, 0.; 1.2, -0.4, 0.1; -0.2, 0.7, 0.5];
    >>> T = [1., 2., -4., 6.; 0., 8., -9., 7.; -3., 7., 8., -2.];
    >>> A2m = pinv(T)*A1m*T;
    >>> B = SimilarityMatrix(A1m, A2m);
    >>> err = norm(A1m*B-B*A2m);
    >>> disp(err);
       1.5605e-15

    For Mathematica:

    >>> A1m = {{0.2, 0.8, 0.},{1.2, -0.4, 0.1},{-0.2, 0.7, 0.5}};
    >>> T = {{1., 2., -4., 6.},{0., 8., -9., 7.},{-3., 7., 8., -2.}};
    >>> A2m = PseudoInverse[T].A1m.T;
    >>> B = SimilarityMatrix[A1m, A2m];
    >>> err = Norm[A1m.B-B.A2m];
    >>> Print[err];
    1.0553656492422553*^-15

    For Python/Numpy:

    >>> A1m = ml.matrix([[0.2, 0.8, 0.],[1.2, -0.4, 0.1],[-0.2, 0.7, 0.5]])
    >>> T = ml.matrix([[1., 2., -4., 6.],[0., 8., -9., 7.],[-3., 7., 8., -2.]])
    >>> A2m = la.pinv(T)*A1m*T
    >>> B = SimilarityMatrix(A1m, A2m)
    >>> err = la.norm(A1m*B-B*A2m)
    >>> print(err)
    9.4758466024e-16

