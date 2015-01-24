butools.dmap.CheckDMRAPRepresentation
=====================================

.. currentmodule:: butools.dmap

.. np:function:: CheckDMRAPRepresentation

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckDMRAPRepresentation(H, prec)`
        * - Mathematica:
          - :code:`r = CheckDMRAPRepresentation[H, prec]`
        * - Python/Numpy:
          - :code:`r = CheckDMRAPRepresentation(H, prec)`

    Checks if the input matrixes define a discrete time MRAP.
    
    All matrices H0...HK must have the same size, the dominant
    eigenvalue of H0 is real and less than 1, and the rowsum of 
    H0+H1+...+HK is 1 (up to the numerical precision).

    Parameters
    ----------
    H : list/cell of matrices, length(K)
        The H0...HK matrices of the DMRAP to check

    Returns
    -------
    r : bool 
        The result of the check

    Examples
    --------
    For Matlab:

    >>> H0=[0.15 0.2 0.18; -0.23 0.17 0.22; 0.19 0.15 0.16];
    >>> H1=[0.01 0.08 0.16; 0.02 0.2 0.07; 0.02 0.15 0.17];
    >>> H2=[0.14 0.07 0.01; 0.19 0.02 0.34; 0.06 0.1 0];
    >>> CheckDMRAPRepresentation({H0,H1,H2})
         1

    For Python/Numpy:
    
    >>> H0=ml.matrix([[0.15, 0.2, 0.18],[-0.23, 0.17, 0.22],[0.19, 0.15, 0.16]])
    >>> H1=ml.matrix([[0.01, 0.08, 0.16],[0.02, 0.2, 0.07],[0.02, 0.15, 0.17]])
    >>> H2=ml.matrix([[0.14, 0.07, 0.01],[0.19, 0.02, 0.34],[0.06, 0.1, 0]])
    >>> print(CheckDMRAPRepresentation((H0,H1,H2)))
    True
    
