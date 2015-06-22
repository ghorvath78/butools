butools.dmap.CheckDMMAPRepresentation
=====================================

.. currentmodule:: butools.dmap

.. np:function:: CheckDMMAPRepresentation

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckDMMAPRepresentation(D, prec)`
        * - Mathematica:
          - :code:`r = CheckDMMAPRepresentation[D, prec]`
        * - Python/Numpy:
          - :code:`r = CheckDMMAPRepresentation(D, prec)`

    Checks if the input matrixes define a discrete time MMAP.
    
    All matrices D0...DK must have the same size, D0 must be a 
    transient probability matrix, D1 has only non-negative 
    elements, and the rowsum of D0+D1+...+DK is 1 (up to the 
    numerical precision).

    Parameters
    ----------
    D : list/cell of matrices, length(K)
        The D0...DK matrices of the DMMAP to check

    Returns
    -------
    r : bool 
        The result of the check

    Examples
    ========
    For Matlab:

    >>> D0 = [0.34, 0, 0; 0.06, 0.05, 0.03; 0.11, 0.13, 0];
    >>> D1 = [0.3, 0, 0; 0.16, 0.18, 0.05; 0.15, 0.04, 0.09];
    >>> D2 = [0, 0.01, 0; 0.1, 0.07, 0.08; 0.13, 0.12, 0.13];
    >>> D3 = [0.35, 0, 0; 0, 0.18, 0.04; 0.06, 0.03, 0.01];
    >>> flag = CheckDMMAPRepresentation({D0,D1,D2,D3});
    >>> disp(flag);
         1

    For Mathematica:

    >>> D0 = {{0.34, 0, 0},{0.06, 0.05, 0.03},{0.11, 0.13, 0}};
    >>> D1 = {{0.3, 0, 0},{0.16, 0.18, 0.05},{0.15, 0.04, 0.09}};
    >>> D2 = {{0, 0.01, 0},{0.1, 0.07, 0.08},{0.13, 0.12, 0.13}};
    >>> D3 = {{0.35, 0, 0},{0, 0.18, 0.04},{0.06, 0.03, 0.01}};
    >>> flag = CheckDMMAPRepresentation[{D0,D1,D2,D3}];
    >>> Print[flag];
    True

    For Python/Numpy:

    >>> D0 = ml.matrix([[0.34, 0, 0],[0.06, 0.05, 0.03],[0.11, 0.13, 0]])
    >>> D1 = ml.matrix([[0.3, 0, 0],[0.16, 0.18, 0.05],[0.15, 0.04, 0.09]])
    >>> D2 = ml.matrix([[0, 0.01, 0],[0.1, 0.07, 0.08],[0.13, 0.12, 0.13]])
    >>> D3 = ml.matrix([[0.35, 0, 0],[0, 0.18, 0.04],[0.06, 0.03, 0.01]])
    >>> flag = CheckDMMAPRepresentation([D0,D1,D2,D3])
    >>> print(flag)
    True

