butools.mc.CheckProbMatrix
==========================

.. currentmodule:: butools.mc

.. np:function:: CheckProbMatrix

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckProbMatrix(P, transient, prec)`
        * - Mathematica:
          - :code:`r = CheckProbMatrix[P, transient, prec]`
        * - Python/Numpy:
          - :code:`r = CheckProbMatrix(P, transient, prec)`
    
    Checks if the matrix is a valid probability matrix: the 
    matrix is a square matrix, the matrix has positive or 
    zero off-diagonal elements, the rowsum of the matrix is 1.
    
    If "transient" is true, it checks if the matrix is a 
    valid transient probability matrix: the matrix is a square
    matrix, the matrix has positive or zero off-diagonal 
    elements, the rowsum of the matrix is less than or equal
    to 1, the maximum absolute eigenvalue is less than 1. 

    Parameters
    ----------
    P : matrix, shape (M,M)
        The matrix to check.
    transient : bool, optional
        If true, the procedure checks if P is a transient 
        probability matrix, otherwise it checks if it is
        a valid probability matrix. The default value is 
        false.
    prec : double, optional
        Entries with absolute value less than prec are 
        considered to be zeros. The default value is 1e-14.
        
    Returns
    -------
    r : bool
        The result of the check.

    Examples
    ========
    For Matlab:

    >>> Q = [0.1, 0.5, 0.4; 0.9, 0.1, 0; 0.3, -0.1, 0.4];
    >>> flag = CheckProbMatrix(Q);
    CheckProbMatrix: the matrix has negative element (precision: 1e-12)!
    >>> disp(flag);
         0
    >>> Q = [0.1, 0.5, 0.4; 0.9, 0.1, 0; 0.3, 0.1, 0.4];
    >>> flag = CheckProbMatrix(Q);
    CheckProbMatrix: The rowsum of the matrix is not 1 (precision: 1e-12)!
    >>> disp(flag);
         0
    >>> Q = [0.1, 0.5, 0.4; 0.9, 0.1, 0; 0.3, 0.3, 0.4];
    >>> flag = CheckProbMatrix(Q);
    >>> disp(flag);
         1
    >>> Q = [0.1, 0.5, 0.4; 0.9, 0.1, 0; 0.3, 0.3, 0.4];
    >>> flag = CheckProbMatrix(Q, true);
    CheckProbMatrix: The real part of the largest eigenvalue of the transient matrix is not less than 1 (precision: 1e-12)!
    >>> disp(flag);
         0
    >>> Q = [0.1, 0.5, 0.4; 0.9, 0.1, 0; 0.3, 0.1, 0.4];
    >>> flag = CheckProbMatrix(Q, true);
    >>> disp(flag);
         1

    For Mathematica:

    >>> Q = {{0.1, 0.5, 0.4},{0.9, 0.1, 0},{0.3, -0.1, 0.4}};
    >>> flag = CheckProbMatrix[Q];
    "CheckProbMatrix: the matrix has negative element (at precision "1.*^-12")!"
    >>> Print[flag];
    False
    >>> Q = {{0.1, 0.5, 0.4},{0.9, 0.1, 0},{0.3, 0.1, 0.4}};
    >>> flag = CheckProbMatrix[Q];
    "CheckProbMatrix: A rowsum of the matrix is not 1 (precision:"1.*^-12")!!"
    >>> Print[flag];
    False
    >>> Q = {{0.1, 0.5, 0.4},{0.9, 0.1, 0},{0.3, 0.3, 0.4}};
    >>> flag = CheckProbMatrix[Q];
    >>> Print[flag];
    True
    >>> Q = {{0.1, 0.5, 0.4},{0.9, 0.1, 0},{0.3, 0.3, 0.4}};
    >>> flag = CheckProbMatrix[Q, True];
    "CheckProbMatrix: The real part of the largest eigenvalue of the transient matrix is not less than 1!"
    >>> Print[flag];
    False
    >>> Q = {{0.1, 0.5, 0.4},{0.9, 0.1, 0},{0.3, 0.1, 0.4}};
    >>> flag = CheckProbMatrix[Q, True];
    >>> Print[flag];
    True

    For Python/Numpy:

    >>> Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, -0.1, 0.4]])
    >>> flag = CheckProbMatrix(Q)
    CheckProbMatrix: the matrix has negative element (precision: 1e-12)!
    >>> print(flag)
    False
    >>> Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, 0.1, 0.4]])
    >>> flag = CheckProbMatrix(Q)
    CheckProbMatrix: The rowsum of the matrix is not 1 (precision: 1e-12)!
    >>> print(flag)
    False
    >>> Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, 0.3, 0.4]])
    >>> flag = CheckProbMatrix(Q)
    >>> print(flag)
    True
    >>> Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, 0.3, 0.4]])
    >>> flag = CheckProbMatrix(Q, True)
    CheckProbMatrix: The real part of the largest eigenvalue of the transient matrix is not less than 1 (precision: 1e-12)!
    >>> print(flag)
    False
    >>> Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, 0.1, 0.4]])
    >>> flag = CheckProbMatrix(Q, True)
    >>> print(flag)
    True

