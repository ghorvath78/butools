butools.map.CheckMRAPRepresentation
===================================

.. currentmodule:: butools.map

.. np:function:: CheckMRAPRepresentation

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckMRAPRepresentation(H, prec)`
        * - Mathematica:
          - :code:`r = CheckMRAPRepresentation[H, prec]`
        * - Python/Numpy:
          - :code:`r = CheckMRAPRepresentation(H, prec)`

    Checks if the input matrixes define a continuous time MRAP.
    
    All matrices H0...HK must have the same size, the dominant
    eigenvalue of H0 is negative and real, and the rowsum of 
    H0+H1+...+HK is 0 (up to the numerical precision).

    Parameters
    ----------
    H : list/cell of matrices, length(K)
        The H0...HK matrices of the MRAP to check

    Returns
    -------
    r : bool 
        The result of the check

    Examples
    ========
    For Matlab:

    >>> x = 0.18;
    >>> H0 = [-5., 0.1+x, 0.9, 1.; 1., -8., 0.9, 0.1; 0.9, 0.1, -4., 1.; 1., 2., 3., -9.];
    >>> H1 = [0.1-x, 0.7, 0.1, 0.1; 0.1, 1., 1.8, 0.1; 0.1, 0.1, 0.1, 0.7; 0.7, 0.1, 0.1, 0.1];
    >>> H2 = [0.1, 0.1, 0.1, 1.7; 1.8, 0.1, 1., 0.1; 0.1, 0.1, 0.7, 0.1; 0.1, 1., 0.1, 0.8];
    >>> flag = CheckMRAPRepresentation({H0, H1, H2});
    >>> disp(flag);
         1

    For Mathematica:

    >>> x = 0.18;
    >>> H0 = {{-5., 0.1+x, 0.9, 1.},{1., -8., 0.9, 0.1},{0.9, 0.1, -4., 1.},{1., 2., 3., -9.}};
    >>> H1 = {{0.1-x, 0.7, 0.1, 0.1},{0.1, 1., 1.8, 0.1},{0.1, 0.1, 0.1, 0.7},{0.7, 0.1, 0.1, 0.1}};
    >>> H2 = {{0.1, 0.1, 0.1, 1.7},{1.8, 0.1, 1., 0.1},{0.1, 0.1, 0.7, 0.1},{0.1, 1., 0.1, 0.8}};
    >>> flag = CheckMRAPRepresentation[{H0, H1, H2}];
    >>> Print[flag];
    True

    For Python/Numpy:

    >>> x = 0.18
    >>> H0 = ml.matrix([[-5., 0.1+x, 0.9, 1.],[1., -8., 0.9, 0.1],[0.9, 0.1, -4., 1.],[1., 2., 3., -9.]])
    >>> H1 = ml.matrix([[0.1-x, 0.7, 0.1, 0.1],[0.1, 1., 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
    >>> H2 = ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1., 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1., 0.1, 0.8]])
    >>> flag = CheckMRAPRepresentation([H0, H1, H2])
    >>> print(flag)
    True

