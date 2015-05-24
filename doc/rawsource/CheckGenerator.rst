butools.mc.CheckGenerator
=========================

.. currentmodule:: butools.mc

.. np:function:: CheckGenerator

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckGenerator(Q, transient, prec)`
        * - Mathematica:
          - :code:`r = CheckGenerator[Q, transient, prec]`
        * - Python/Numpy:
          - :code:`r = CheckGenerator(Q, transient, prec)`
    
    Checks if the matrix is a valid generator matrix: the 
    matrix is a square matrix, the matrix has positive or 
    zero off-diagonal elements, the diagonal of the matrix 
    is negative, the rowsum of the matrix is 0.
 
    If the "transient" parameter is set to false, it checks 
    if the real part of the maximum absolute eigenvalue is 
    less than zero and the rowsum is equal or less than 0. 

    Parameters
    ----------
    Q : matrix, shape (M,M)
        The generator to check.
    transient : bool, optional
        If true, the procedure checks if Q is a transient 
        generator, otherwise it checks if it is a valid 
        generator. The default value is false.
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
    
    >>> q=[-0.9 0.2 0.4; 0 -0.9 0.9; 0 0.6 -0.6]
    >>> CheckGenerator(q, true)
    1
    >>> CheckGenerator(q)
    0

    For Mathematica:
    
    >>> q={{-0.9, 0.2, 0.4}, {0, -0.9, 0.9}, {0, 0.6, -0.6}};
    >>> CheckGenerator[q, True]
    True
    >>> CheckGenerator[q]
    False    
    
    For Python/Numpy:
    
    >>> Q=[[-0.9, 0.2, 0.4], [0, -0.9, 0.9], [0, 0.6, -0.6]]
    >>> CheckGenerator(Q, True)
    True
    >>> CheckGenerator(Q)
    False    
    
