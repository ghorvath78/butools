butools.dmap.CheckDMAPRepresentation
====================================

.. currentmodule:: butools.dmap

.. np:function:: CheckDMAPRepresentation

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckDMAPRepresentation(D0, D1, prec)`
        * - Mathematica:
          - :code:`r = CheckDMAPRepresentation[D0, D1, prec]`
        * - Python/Numpy:
          - :code:`r = CheckDMAPRepresentation(D0, D1, prec)`

    Checks if the input matrixes define a discrete time MAP.

    Matrices D0 and D1 must have the same size, D0 must be a 
    transient probability matrix, D1 has only non-negative
    elements, and the rowsum of D0+D1 is 1 (up to the numerical
    precision).

    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the DMAP to check
    D1 : matrix, shape (M,M)
        The D1 matrix of the DMAP to check
    prec : double, optional
        Numerical precision, the default value is 1e-14

    Returns
    -------
    r : bool 
        The result of the check

    Examples
    --------
    For Matlab:
    
    >>> D0=[0 0.02 0; 0 0.17 0.2; 0.16 0.17 0.02]
    >>> D1=[0 0.88 0.1 0; 0.18 0.07 0.14 0.1; 0.13 0.15 0.15 0.04; 0.31 0.18 0.12 0.27]
    >>> CheckDMAPRepresentation(D0,D1)
        CheckDMAPRepresentation: D0 and D1 have different sizes!
             0
    >>> D0=[0 0.02 0; 0 0.17 0.2; 0.16 0.17 0.02]
    >>> D1=[0 0.88 0.1; 0.18 0.07 0.14; 0.13 0.15 0.15]
    >>> CheckDMAPRepresentation(D0,D1)
        CheckDMAPRepresentation: A rowsum of matrix0+matrix1 is not 1 (at precision 1e-14)!
             0
    >>> D0=[0 0.02 0 0; 0 0.17 0.2 0.14; 0.16 0.17 0.02 0.18; 0 0 0 0.12];
    >>> D1=[0 0.88 0.1 0; 0.18 0.07 0.14 0.1; 0.13 0.15 0.15 0.04; 0.31 0.18 0.12 0.27];
    >>> CheckDMAPRepresentation(D0,D1)
             1

