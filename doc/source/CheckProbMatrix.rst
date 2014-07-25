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
    --------
    For Matlab:
    
    >>> q=[0.1 0.5 0.4; 0.9 0.1 0; 0.3 0.1 0.4]
    >>> CheckProbMatrix(q, true)
    1
    >>> CheckProbMatrix(q)
    0

    For Mathematica:
    
    >>> q={{0.1, 0.5, 0.4}, {0.9, 0.1, 0}, {0.3, 0.1, 0.4}}
    >>> CheckProbMatrix[q, True]
    >>> CheckProbMatrix[q]
    
    For Python/Numpy:
    
    >>> q= [[0.1, 0.5, 0.4], [0.9, 0.1, 0], [0.3, 0.1, 0.4]]
    >>> CheckProbMatrix(q, True)
    >>> CheckProbMatrix(q)

