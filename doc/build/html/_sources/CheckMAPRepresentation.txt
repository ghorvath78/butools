butools.map.CheckMAPRepresentation
==================================

.. currentmodule:: butools.map

.. np:function:: CheckMAPRepresentation

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckMAPRepresentation(D0, D1, prec)`
        * - Mathematica:
          - :code:`r = CheckMAPRepresentation[D0, D1, prec]`
        * - Python/Numpy:
          - :code:`r = CheckMAPRepresentation(D0, D1, prec)`

    Checks if the input matrixes define a continuous time MAP.
    
    Matrices D0 and D1 must have the same size, D0 must be a 
    transient generator matrix, D1 has only non-negative 
    elements, and the rowsum of D0+D1 is 0 (up to the numerical
    precision).

    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the MAP to check
    D1 : matrix, shape (M,M)
        The D1 matrix of the MAP to check
    prec : double, optional
        Numerical precision, the default value is 1e-14

    Returns
    -------
    r : bool 
        The result of the check

    Examples
    ========
    For Matlab:

    >>> D0 = [-1., 0, 1.; 0, -2., 0; 1., 0, -3.];
    >>> D1 = [-1., 0, 1., 0; 0, -2., 0, 1.; 1., 0, -3., 0; 1., 2., 2., 1.];
    >>> flag = CheckMAPRepresentation(D0,D1);
    CheckMAPRepresentation: D0 and D1 have different sizes!
    >>> disp(flag);
         0
    >>> D0 = [-1., 0, 1.; 0, -2., 0; 1., 0, -3.];
    >>> D1 = [1., 0, 1.; 0, 2., 0; 1., 0, 3.];
    >>> flag = CheckMAPRepresentation(D0,D1);
    CheckMAPRepresentation: The rowsum of D0+D1 is not 0 (precision: 1e-10)!
    >>> disp(flag);
         0
    >>> D0 = [-3., 0, 1.; 0, -2., 0; 1., 0, -5.];
    >>> D1 = [1., 0, 1.; 0, 2., 0; 1., 0, 3.];
    >>> flag = CheckMAPRepresentation(D0,D1);
    >>> disp(flag);
         1

    For Mathematica:

    >>> D0 = {{-1., 0, 1.},{0, -2., 0},{1., 0, -3.}};
    >>> D1 = {{-1., 0, 1., 0},{0, -2., 0, 1.},{1., 0, -3., 0},{1., 2., 2., 1.}};
    >>> flag = CheckMAPRepresentation[D0,D1];
    "CheckMAPRepresentation: D0 and D1 have different sizes!"
    >>> Print[flag];
    False
    >>> D0 = {{-1., 0, 1.},{0, -2., 0},{1., 0, -3.}};
    >>> D1 = {{1., 0, 1.},{0, 2., 0},{1., 0, 3.}};
    >>> flag = CheckMAPRepresentation[D0,D1];
    "CheckMAPRepresentation: A rowsum of D0+D1 is not 0 (at precision "1.*^-12")!"
    >>> Print[flag];
    False
    >>> D0 = {{-3., 0, 1.},{0, -2., 0},{1., 0, -5.}};
    >>> D1 = {{1., 0, 1.},{0, 2., 0},{1., 0, 3.}};
    >>> flag = CheckMAPRepresentation[D0,D1];
    >>> Print[flag];
    True

    For Python/Numpy:

    >>> D0 = ml.matrix([[-1., 0, 1.],[0, -2., 0],[1., 0, -3.]])
    >>> D1 = ml.matrix([[-1., 0, 1., 0],[0, -2., 0, 1.],[1., 0, -3., 0],[1., 2., 2., 1.]])
    >>> flag = CheckMAPRepresentation(D0,D1)
    CheckMAPRepresentation: D0 and D1 have different sizes!
    >>> print(flag)
    False
    >>> D0 = ml.matrix([[-1., 0, 1.],[0, -2., 0],[1., 0, -3.]])
    >>> D1 = ml.matrix([[1., 0, 1.],[0, 2., 0],[1., 0, 3.]])
    >>> flag = CheckMAPRepresentation(D0,D1)
    CheckMAPRepresentation: The rowsum of D0+D1 is not 0!
    >>> print(flag)
    False
    >>> D0 = ml.matrix([[-3., 0, 1.],[0, -2., 0],[1., 0, -5.]])
    >>> D1 = ml.matrix([[1., 0, 1.],[0, 2., 0],[1., 0, 3.]])
    >>> flag = CheckMAPRepresentation(D0,D1)
    >>> print(flag)
    True

