butools.mc.CheckProbVector
==========================

.. currentmodule:: butools.mc

.. np:function:: CheckProbVector

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckProbVector(pi, sub, prec)`
        * - Mathematica:
          - :code:`r = CheckProbVector[pi, sub, prec]`
        * - Python/Numpy:
          - :code:`r = CheckProbVector(pi, sub, prec)`
    
    Checks if the vector is a valid probability vector: the 
    vector has only non-negative elements, the sum of the 
    vector elements is 1.
 
    If parameter "sub" is set to true, it checks if the 
    vector is a valid substochastic vector: the vector has 
    only non-negative elements, the sum of the elements are
    less than 1.
    
    Parameters
    ----------
    pi : vector, shape (1, M) or (M, 1)
        The matrix to check.
    sub : bool, optional
        If false, the procedure checks for stochastic, if 
        true, it checks for sub-stochastic property. The 
        default value is false.
    prec : double, optional
        Numerical precision. Entries with absolute value 
        less than prec are considered to be zeros. The 
        default value is 1e-14.
        
    Returns
    -------
    r : bool
        The result of the check.

    Examples
    ========
    For Matlab:

    >>> Q = [1.1, -0.1];
    >>> flag = CheckProbVector(Q);
    CheckProbVector: The vector has negative element (precision: 1e-12)!
    >>> disp(flag);
         0
    >>> Q = [1.1, 0.1];
    >>> flag = CheckProbVector(Q);
    CheckProbVector: The sum of the vector is not 1 (precision: 1e-12)!
    >>> disp(flag);
         0
    >>> Q = [1, 0];
    >>> flag = CheckProbVector(Q);
    >>> disp(flag);
         1
    >>> Q = [0.9, -0.1];
    >>> flag = CheckProbVector(Q,true);
    CheckProbVector: The vector has negative element (precision: 1e-12)!
    >>> disp(flag);
         0
    >>> Q = [0.9, 0.1];
    >>> flag = CheckProbVector(Q,true);
    >>> disp(flag);
         1
    >>> Q = [0.8, 0.1];
    >>> flag = CheckProbVector(Q,true);
    >>> disp(flag);
         1

    For Mathematica:

    >>> Q = {1.1, -0.1};
    >>> flag = CheckProbVector[Q];
    "CheckProbVector: The vector has negative element!"
    >>> Print[flag];
    False
    >>> Q = {1.1, 0.1};
    >>> flag = CheckProbVector[Q];
    "CheckProbVector: The sum of the vector is not 1 (precision:"1.*^-12")!"
    >>> Print[flag];
    False
    >>> Q = {1, 0};
    >>> flag = CheckProbVector[Q];
    >>> Print[flag];
    True
    >>> Q = {0.9, -0.1};
    >>> flag = CheckProbVector[Q,True];
    "CheckProbVector: The vector has negative element!"
    >>> Print[flag];
    False
    >>> Q = {0.9, 0.1};
    >>> flag = CheckProbVector[Q,True];
    >>> Print[flag];
    True
    >>> Q = {0.8, 0.1};
    >>> flag = CheckProbVector[Q,True];
    >>> Print[flag];
    True

    For Python/Numpy:

    >>> Q = ml.matrix([[1.1, -0.1]])
    >>> flag = CheckProbVector(Q)
    CheckProbVector: The vector has negative element (precision: 1e-12)!
    >>> print(flag)
    False
    >>> Q = ml.matrix([[1.1, 0.1]])
    >>> flag = CheckProbVector(Q)
    CheckProbVector: The sum of the vector is not 1 (precision: 1e-12)!
    >>> print(flag)
    False
    >>> Q = ml.matrix([[1, 0]])
    >>> flag = CheckProbVector(Q)
    >>> print(flag)
    True
    >>> Q = ml.matrix([[0.9, -0.1]])
    >>> flag = CheckProbVector(Q,True)
    CheckProbVector: The vector has negative element (precision: 1e-12)!
    >>> print(flag)
    False
    >>> Q = ml.matrix([[0.9, 0.1]])
    >>> flag = CheckProbVector(Q,True)
    >>> print(flag)
    True
    >>> Q = ml.matrix([[0.8, 0.1]])
    >>> flag = CheckProbVector(Q,True)
    >>> print(flag)
    True

