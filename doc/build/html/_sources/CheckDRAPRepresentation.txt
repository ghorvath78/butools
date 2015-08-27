butools.dmap.CheckDRAPRepresentation
====================================

.. currentmodule:: butools.dmap

.. np:function:: CheckDRAPRepresentation

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckDRAPRepresentation(H, prec)`
        * - Mathematica:
          - :code:`r = CheckDRAPRepresentation[H, prec]`
        * - Python/Numpy:
          - :code:`r = CheckDRAPRepresentation(H, prec)`

    Checks if the input matrixes define a discrete time RAP.
    
    Matrices H0 and H1 must have the same size, the dominant
    eigenvalue of H0 is real and less than 1, and the rowsum of 
    H0+H1 is 1 (up to the numerical precision).

    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the DRAP to check
    H1 : matrix, shape (M,M)
        The H1 matrix of the DRAP to check
    prec : double, optional
        Numerical precision, the default value is 1e-14

    Returns
    -------
    r : bool 
        The result of the check

    Examples
    ========
    For Matlab:

    >>> H0 = [0, 0, 0.13; 0, 0.6, 0.18; 0.31, 0.26, 0.02; 0.2, 0, 0];
    >>> H1 = [0, 1., -0.13; 0, 0.18, 0.04; 0.03, 0.09, 0.29; 0, 0.8, 0];
    >>> flag = CheckDRAPRepresentation(H0, H1);
    CheckDRAPRepresentation: D0 and D1 have different sizes!
    >>> disp(flag);
         0
    >>> H0 = [0.2, 0, 0.13; 0, 0.6, 0.18; 0.31, 0.26, 0.02];
    >>> H1 = [0, 1., -0.13; 0, 0.18, 0.04; 0.03, 0.09, 0.29];
    >>> flag = CheckDRAPRepresentation(H0, H1);
    CheckDRAPRepresentation: A rowsum of D0+D1 is not 1 (at precision 1e-12)!
    >>> disp(flag);
         0
    >>> H0 = [-1., 0, 0; 0, -2., 2.; 0, 3., -3.];
    >>> H1 = [0, 0, 1.; 0, -1., 1.; 1., 0, -1.];
    >>> flag = CheckDRAPRepresentation(H0, H1);
    CheckDRAPRepresentation: A rowsum of D0+D1 is not 1 (at precision 1e-12)!
    >>> disp(flag);
         0
    >>> H0 = [0, 0, 15.; 0, 0.6, 0.18; 0.31, 0.26, 0.02];
    >>> H1 = [0, 1., -15.; 0, 0.18, 0.04; 0.03, 0.09, 0.29];
    >>> flag = CheckDRAPRepresentation(H0, H1);
    CheckDRAPRepresentation: The dominant eigenvalue of D0 is greater than 1!
    >>> disp(flag);
         0
    >>> H0 = [0, 0.5, 0.1; 0, -1.4, 3.1; 0.67, 0, 0.4];
    >>> H1 = [0, 0.4, 0; 0, -0.2, -0.5; 0.3, -0.7, 0.33];
    >>> flag = CheckDRAPRepresentation(H0, H1);
    CheckDRAPRepresentation: The dominant eigenvalue of the D0 is complex!
    >>> disp(flag);
         0
    >>> H0 = [0, 0, 0.13; 0, 0.6, 0.18; 0.31, 0.26, 0.02];
    >>> H1 = [0, 1., -0.13; 0, 0.18, 0.04; 0.03, 0.09, 0.29];
    >>> flag = CheckDRAPRepresentation(H0, H1);
    >>> disp(flag);
         1

    For Mathematica:

    >>> H0 = {{0, 0, 0.13},{0, 0.6, 0.18},{0.31, 0.26, 0.02},{0.2, 0, 0}};
    >>> H1 = {{0, 1., -0.13},{0, 0.18, 0.04},{0.03, 0.09, 0.29},{0, 0.8, 0}};
    >>> flag = CheckDRAPRepresentation[H0, H1];
    "CheckDRAPRepresentation: D0 is not a quadratic matrix!"
    >>> Print[flag];
    False
    >>> H0 = {{0.2, 0, 0.13},{0, 0.6, 0.18},{0.31, 0.26, 0.02}};
    >>> H1 = {{0, 1., -0.13},{0, 0.18, 0.04},{0.03, 0.09, 0.29}};
    >>> flag = CheckDRAPRepresentation[H0, H1];
    "CheckDRAPRepresentation: A rowsum of D0+D1 is not 1! (precision:"1.*^-12")"
    >>> Print[flag];
    False
    >>> H0 = {{-1., 0, 0},{0, -2., 2.},{0, 3., -3.}};
    >>> H1 = {{0, 0, 1.},{0, -1., 1.},{1., 0, -1.}};
    >>> flag = CheckDRAPRepresentation[H0, H1];
    "CheckDRAPRepresentation: A rowsum of D0+D1 is not 1! (precision:"1.*^-12")"
    >>> Print[flag];
    False
    >>> H0 = {{0, 0, 15.},{0, 0.6, 0.18},{0.31, 0.26, 0.02}};
    >>> H1 = {{0, 1., -15.},{0, 0.18, 0.04},{0.03, 0.09, 0.29}};
    >>> flag = CheckDRAPRepresentation[H0, H1];
    "CheckDRAPRepresentation: The dominant eigenvalue of D0 is greater than 1!"
    >>> Print[flag];
    False
    >>> H0 = {{0, 0.5, 0.1},{0, -1.4, 3.1},{0.67, 0, 0.4}};
    >>> H1 = {{0, 0.4, 0},{0, -0.2, -0.5},{0.3, -0.7, 0.33}};
    >>> flag = CheckDRAPRepresentation[H0, H1];
    "CheckDRAPRepresentation: The dominant eigenvalue of D0 is complex!"
    >>> Print[flag];
    False
    >>> H0 = {{0, 0, 0.13},{0, 0.6, 0.18},{0.31, 0.26, 0.02}};
    >>> H1 = {{0, 1., -0.13},{0, 0.18, 0.04},{0.03, 0.09, 0.29}};
    >>> flag = CheckDRAPRepresentation[H0, H1];
    >>> Print[flag];
    True

    For Python/Numpy:

    >>> H0 = ml.matrix([[0, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02],[0.2, 0, 0]])
    >>> H1 = ml.matrix([[0, 1., -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29],[0, 0.8, 0]])
    >>> flag = CheckDRAPRepresentation(H0, H1)
    CheckDRAPRepresentation: D0 is not a quadratic matrix!
    >>> print(flag)
    False
    >>> H0 = ml.matrix([[0.2, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    >>> H1 = ml.matrix([[0, 1., -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])
    >>> flag = CheckDRAPRepresentation(H0, H1)
    CheckDRAPRepresentation: A rowsum of D0+D1 is not 1!
    >>> print(flag)
    False
    >>> H0 = ml.matrix([[-1., 0, 0],[0, -2., 2.],[0, 3., -3.]])
    >>> H1 = ml.matrix([[0, 0, 1.],[0, -1., 1.],[1., 0, -1.]])
    >>> flag = CheckDRAPRepresentation(H0, H1)
    CheckDRAPRepresentation: A rowsum of D0+D1 is not 1!
    >>> print(flag)
    False
    >>> H0 = ml.matrix([[0, 0, 15.],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    >>> H1 = ml.matrix([[0, 1., -15.],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])
    >>> flag = CheckDRAPRepresentation(H0, H1)
    CheckDRAPRepresentation: The largest eigenvalue of matrix D0 is greater than 1!
    >>> print(flag)
    False
    >>> H0 = ml.matrix([[0, 0.5, 0.1],[0, -1.4, 3.1],[0.67, 0, 0.4]])
    >>> H1 = ml.matrix([[0, 0.4, 0],[0, -0.2, -0.5],[0.3, -0.7, 0.33]])
    >>> flag = CheckDRAPRepresentation(H0, H1)
    CheckDRAPRepresentation: The largest eigenvalue of matrix D0 is complex!
    >>> print(flag)
    False
    >>> H0 = ml.matrix([[0, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    >>> H1 = ml.matrix([[0, 1., -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])
    >>> flag = CheckDRAPRepresentation(H0, H1)
    >>> print(flag)
    True

