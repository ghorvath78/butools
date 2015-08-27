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
    ========
    For Matlab:

    >>> Q = [-0.9, 0.2, 0.4; 0, 0.9, 0.9; 0, 0.6, -0.6];
    >>> flag = CheckGenerator(Q, true);
    CheckGenerator: The diagonal of the generator is not negative (precision: 1e-12)!
    >>> disp(flag);
         0
    >>> Q = [-0.9, 0.5, 0.4; 0.9, -0.9, 0; 0.3, 0.3, -0.6];
    >>> flag = CheckGenerator(Q, true);
    >>> disp(flag);
         1
    >>> Q = [-0.9, 0.2, 0.4; 0.9, -0.9, 0; 0.3, 0.3, -0.6];
    >>> flag = CheckGenerator(Q, true);
    >>> disp(flag);
         1
    >>> Q = [-0.9, 0.5, 0.4; 0.9, -1.1, 0; 0.3, 0.3, -0.6];
    >>> flag = CheckGenerator(Q);
    CheckGenerator: The rowsum of the generator is not 0 (precision: 1e-12)!
    >>> disp(flag);
         0
    >>> Q = [-0.9, 0.5, 0.4; 0.9, -0.9, 0; 0.3, 0.3, -0.6];
    >>> flag = CheckGenerator(Q);
    >>> disp(flag);
         1

    For Mathematica:

    >>> Q = {{-0.9, 0.2, 0.4},{0, 0.9, 0.9},{0, 0.6, -0.6}};
    >>> flag = CheckGenerator[Q, True];
    "CheckGenerator: The diagonal of the generator is not negative (at precision "1.*^-12")!"
    >>> Print[flag];
    False
    >>> Q = {{-0.9, 0.5, 0.4},{0.9, -0.9, 0},{0.3, 0.3, -0.6}};
    >>> flag = CheckGenerator[Q, True];
    >>> Print[flag];
    True
    >>> Q = {{-0.9, 0.2, 0.4},{0.9, -0.9, 0},{0.3, 0.3, -0.6}};
    >>> flag = CheckGenerator[Q, True];
    >>> Print[flag];
    True
    >>> Q = {{-0.9, 0.5, 0.4},{0.9, -1.1, 0},{0.3, 0.3, -0.6}};
    >>> flag = CheckGenerator[Q];
    "CheckGenerator: A rowsum of the generator is not 0 (precision:"1.*^-12")!!"
    >>> Print[flag];
    False
    >>> Q = {{-0.9, 0.5, 0.4},{0.9, -0.9, 0},{0.3, 0.3, -0.6}};
    >>> flag = CheckGenerator[Q];
    >>> Print[flag];
    True

    For Python/Numpy:

    >>> Q = ml.matrix([[-0.9, 0.2, 0.4],[0, 0.9, 0.9],[0, 0.6, -0.6]])
    >>> flag = CheckGenerator(Q, True)
    CheckGenerator: The diagonal of the generator is not negative (precision: 1e-12)!
    >>> print(flag)
    False
    >>> Q = ml.matrix([[-0.9, 0.5, 0.4],[0.9, -0.9, 0],[0.3, 0.3, -0.6]])
    >>> flag = CheckGenerator(Q, True)
    >>> print(flag)
    True
    >>> Q = ml.matrix([[-0.9, 0.2, 0.4],[0.9, -0.9, 0],[0.3, 0.3, -0.6]])
    >>> flag = CheckGenerator(Q, True)
    >>> print(flag)
    True
    >>> Q = ml.matrix([[-0.9, 0.5, 0.4],[0.9, -1.1, 0],[0.3, 0.3, -0.6]])
    >>> flag = CheckGenerator(Q)
    CheckGenerator: The rowsum of the generator is not 0 (precision: 1e-12)!
    >>> print(flag)
    False
    >>> Q = ml.matrix([[-0.9, 0.5, 0.4],[0.9, -0.9, 0],[0.3, 0.3, -0.6]])
    >>> flag = CheckGenerator(Q)
    >>> print(flag)
    True

