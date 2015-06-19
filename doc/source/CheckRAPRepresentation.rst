butools.map.CheckRAPRepresentation
==================================

.. currentmodule:: butools.map

.. np:function:: CheckRAPRepresentation

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckRAPRepresentation(H0, H1, prec)`
        * - Mathematica:
          - :code:`r = CheckRAPRepresentation[H0, H1, prec]`
        * - Python/Numpy:
          - :code:`r = CheckRAPRepresentation(H0, H1, prec)`

    Checks if the input matrixes define a continuous time RAP.
    
    Matrices H0 and H1 must have the same size, the dominant
    eigenvalue of H0 is negative and real, and the rowsum of 
    H0+H1 is 0 (up to the numerical precision).

    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the RAP to check
    H1 : matrix, shape (M,M)
        The H1 matrix of the RAP to check
    prec : double, optional
        Numerical precision, the default value is 1e-14

    Returns
    -------
    r : bool 
        The result of the check

    Examples
    ========
    For Matlab:

    >>> H0 = [-1., 0, 1.; 0, -2., 0; 1., 0, -3.; 1., 2., 2.];
    >>> H1 = [-1., 0, 1.; 0, -2., 0; 1., 0, -3.; 1., 2., 2.];
    >>> flag = CheckRAPRepresentation(H0,H1);
    CheckRAPRepresentation: D0 is not a quadratic matrix!
    >>> disp(flag);
         0
    >>> H0 = [-1., 0, 2.; 0, 2., 0; 1., 0, -3.];
    >>> H1 = [-1., 0, 1.; 0, -2., 0; 1., 0, -3.];
    >>> flag = CheckRAPRepresentation(H0,H1);
    CheckRAPRepresentation: A rowsum of D0+D1 is not 0!(precision: 1e-11)
    >>> disp(flag);
         0
    >>> H0 = [-1., 0, 0; 0, -2., 2.; 0, 3., -3.];
    >>> H1 = [0, 0, 1.; 0, -1., 1.; 1., 0, -1.];
    >>> flag = CheckRAPRepresentation(H0,H1);
    CheckRAPRepresentation: there is an eigenvalue of D0 with non-negative real part (at precision 1e-11)
    >>> disp(flag);
         0
    >>> H0 = [-2., 0, 0; 0, -1., 1.; 0, -1., -1.];
    >>> H1 = [1., 0, 1.; 0, 1., -1.; 1., 0, 1.];
    >>> flag = CheckRAPRepresentation(H0,H1);
    CheckRAPRepresentation: The dominant eigenvalue of D0 is not real!
    >>> disp(flag);
         0
    >>> H0 = [-1., 0, 0; 0, -2., 1.; 0, -1., -2.];
    >>> H1 = [1., 0, 0; 0, 1., 0; 1., 1., 1.];
    >>> flag = CheckRAPRepresentation(H0,H1);
    >>> disp(flag);
         1

    For Mathematica:

    >>> H0 = {{-1., 0, 1.},{0, -2., 0},{1., 0, -3.},{1., 2., 2.}};
    >>> H1 = {{-1., 0, 1.},{0, -2., 0},{1., 0, -3.},{1., 2., 2.}};
    >>> flag = CheckRAPRepresentation[H0,H1];
    "CheckRAPRepresentation: D0 is not a quadratic matrix!"
    >>> Print[flag];
    False
    >>> H0 = {{-1., 0, 2.},{0, 2., 0},{1., 0, -3.}};
    >>> H1 = {{-1., 0, 1.},{0, -2., 0},{1., 0, -3.}};
    >>> flag = CheckRAPRepresentation[H0,H1];
    "CheckRAPRepresentation: A rowsum of D0+D1 is not 0! (precision:"1.*^-12")"
    >>> Print[flag];
    False
    >>> H0 = {{-1., 0, 0},{0, -2., 2.},{0, 3., -3.}};
    >>> H1 = {{0, 0, 1.},{0, -1., 1.},{1., 0, -1.}};
    >>> flag = CheckRAPRepresentation[H0,H1];
    "CheckRAPRepresentation: there is an eigenvalue of D0 with non-negative real part (at precision "1.*^-12")!"
    >>> Print[flag];
    False
    >>> H0 = {{-2., 0, 0},{0, -1., 1.},{0, -1., -1.}};
    >>> H1 = {{1., 0, 1.},{0, 1., -1.},{1., 0, 1.}};
    >>> flag = CheckRAPRepresentation[H0,H1];
    "CheckRAPRepresentation: The dominant eigenvalue of D0 is not real!"
    >>> Print[flag];
    False
    >>> H0 = {{-1., 0, 0},{0, -2., 1.},{0, -1., -2.}};
    >>> H1 = {{1., 0, 0},{0, 1., 0},{1., 1., 1.}};
    >>> flag = CheckRAPRepresentation[H0,H1];
    >>> Print[flag];
    True

    For Python/Numpy:

    >>> H0 = ml.matrix([[-1., 0, 1.],[0, -2., 0],[1., 0, -3.],[1., 2., 2.]])
    >>> H1 = ml.matrix([[-1., 0, 1.],[0, -2., 0],[1., 0, -3.],[1., 2., 2.]])
    >>> flag = CheckRAPRepresentation(H0,H1)
    CheckRAPRepresentation: D0 is not a quadratic matrix!
    >>> print(flag)
    False
    >>> H0 = ml.matrix([[-1., 0, 2.],[0, 2., 0],[1., 0, -3.]])
    >>> H1 = ml.matrix([[-1., 0, 1.],[0, -2., 0],[1., 0, -3.]])
    >>> flag = CheckRAPRepresentation(H0,H1)
    CheckRAPRepresentation: A rowsum of D0+D1 is not 0!
    >>> print(flag)
    False
    >>> H0 = ml.matrix([[-1., 0, 0],[0, -2., 2.],[0, 3., -3.]])
    >>> H1 = ml.matrix([[0, 0, 1.],[0, -1., 1.],[1., 0, -1.]])
    >>> flag = CheckRAPRepresentation(H0,H1)
    CheckRAPRepresentation: there is an eigenvalue of D0 with non-negative real part
    >>> print(flag)
    False
    >>> H0 = ml.matrix([[-2., 0, 0],[0, -1., 1.],[0, -1., -1.]])
    >>> H1 = ml.matrix([[1., 0, 1.],[0, 1., -1.],[1., 0, 1.]])
    >>> flag = CheckRAPRepresentation(H0,H1)
    CheckRAPRepresentation: The dominant eigenvalue of D0 is not real!
    >>> print(flag)
    False
    >>> H0 = ml.matrix([[-1., 0, 0],[0, -2., 1.],[0, -1., -2.]])
    >>> H1 = ml.matrix([[1., 0, 0],[0, 1., 0],[1., 1., 1.]])
    >>> flag = CheckRAPRepresentation(H0,H1)
    >>> print(flag)
    True

