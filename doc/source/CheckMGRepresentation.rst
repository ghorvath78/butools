butools.dph.CheckMGRepresentation
=================================

.. currentmodule:: butools.dph

.. np:function:: CheckMGRepresentation

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckMGRepresentation(alpha, A, prec)`
        * - Mathematica:
          - :code:`r = CheckMGRepresentation[alpha, A, prec]`
        * - Python/Numpy:
          - :code:`r = CheckMGRepresentation(alpha, A, prec)`

    Checks if the given vector and matrix define a valid matrix-
    geometric representation.

    Parameters
    ----------
    alpha : matrix, shape (1,M)
        Initial vector of the matrix-geometric distribution 
        to check
    A : matrix, shape (M,M)
        Matrix parameter of the matrix-geometric distribution
        to check
    prec : double, optional
        Numerical precision. The default value is 1e-14.

    Returns
    -------
    r : bool
        True, if the matrix is a square matrix, the vector and 
        the matrix have the same size, the dominant eigenvalue
        is positive, less than 1 and real. 

    Notes
    -----
    This procedure does not check the positivity of the density!
    The discrete counterpart of 'CheckMEPositiveDensity' does
    not exist yet (research is needed).

    Examples
    ========
    For Matlab:

    >>> a = [-0.6,0.3,1.3];
    >>> A = [0.25, 0.2, -0.15; 0.3, 0.1, 0.25; 0, 0.2, 0.47];
    >>> flag = CheckMGRepresentation(a, A);
    >>> disp(flag);
         1
    >>> a = [-0.6,0.3,1.3];
    >>> A = [0.35, 0.2, -0.25; 0.3, 0.1, 0.25; 0, 0.2, 0.47];
    >>> flag = CheckMGRepresentation(a, A);
    CheckMGRepresentation: The largest eigenvalue of the matrix is complex!
    >>> disp(flag);
         0

    For Mathematica:

    
    For Python/Numpy:

    >>> a = ml.matrix([[-0.6,0.3,1.3]])
    >>> A = ml.matrix([[0.25, 0.2, -0.15],[0.3, 0.1, 0.25],[0, 0.2, 0.47]])
    >>> flag = CheckMGRepresentation(a, A)
    >>> print(flag)
    True
    >>> a = ml.matrix([[-0.6,0.3,1.3]])
    >>> A = ml.matrix([[0.35, 0.2, -0.25],[0.3, 0.1, 0.25],[0, 0.2, 0.47]])
    >>> flag = CheckMGRepresentation(a, A)
    CheckMGRepresentation: The largest eigenvalue of the matrix is complex!
    >>> print(flag)
    False

