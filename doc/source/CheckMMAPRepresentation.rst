butools.map.CheckMMAPRepresentation
===================================

.. currentmodule:: butools.map

.. np:function:: CheckMMAPRepresentation

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckMMAPRepresentation(D, prec)`
        * - Mathematica:
          - :code:`r = CheckMMAPRepresentation[D, prec]`
        * - Python/Numpy:
          - :code:`r = CheckMMAPRepresentation(D, prec)`

    Checks if the input matrixes define a continuous time MMAP.
    
    All matrices D0...DK must have the same size, D0 must be a 
    transient generator matrix, D1 has only non-negative 
    elements, and the rowsum of D0+D1+...+DK is 0 (up to the 
    numerical precision).

    Parameters
    ----------
    D : list/cell of matrices, length(K)
        The D0...DK matrices of the MMAP to check

    Returns
    -------
    r : bool 
        The result of the check

    Examples
    ========
    For Matlab:

    >>> D0 = [-1.05, 0.03, 0.07; 0.19, -1.63, 0.06; 0, 0.2, -1.03];
    >>> D1 = [0.16, 0.11, 0; 0.1, 0.16, 0; 0.27, 0, 0.19];
    >>> D2 = [0.01, 0.09, 0.13; 0.26, 0.21, 0.05; 0, 0.16, 0.07];
    >>> D3 = [0.19, 0.06, 0.2; 0.17, 0.16, 0.27; 0, 0, 0.14];
    >>> flag = CheckMMAPRepresentation({D0, D1, D2, D3});
    >>> disp(flag);
         1

    For Mathematica:

    
    For Python/Numpy:

    >>> D0 = ml.matrix([[-1.05, 0.03, 0.07],[0.19, -1.63, 0.06],[0, 0.2, -1.03]])
    >>> D1 = ml.matrix([[0.16, 0.11, 0],[0.1, 0.16, 0],[0.27, 0, 0.19]])
    >>> D2 = ml.matrix([[0.01, 0.09, 0.13],[0.26, 0.21, 0.05],[0, 0.16, 0.07]])
    >>> D3 = ml.matrix([[0.19, 0.06, 0.2],[0.17, 0.16, 0.27],[0, 0, 0.14]])
    >>> flag = CheckMMAPRepresentation([D0, D1, D2, D3])
    >>> print(flag)
    True

