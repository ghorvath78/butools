butools.ph.CheckPHRepresentation
================================

.. currentmodule:: butools.ph

.. np:function:: CheckPHRepresentation

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckPHRepresentation(alpha, A, prec)`
        * - Mathematica:
          - :code:`r = CheckPHRepresentation[alpha, A, prec]`
        * - Python/Numpy:
          - :code:`r = CheckPHRepresentation(alpha, A, prec)`

    Checks if the given vector and matrix define a valid phase-
    type representation.

    Parameters
    ----------
    alpha : matrix, shape (1,M)
        Initial vector of the phase-type distribution to check
    A : matrix, shape (M,M)
        Transient generator of the phase-type distribution to
        check
    prec : double, optional
        Numerical precision. The default value is 1e-14.

    Returns
    -------
    r : bool
        True, if vector alpha is a probability vector and matrix
        A is a transient generator, and they have the same size.

    Examples
    --------
    For Matlab:
    
    >>> a=[0.2];
    >>> A=[-1 1; 1 -2];
    >>> CheckPHRepresentation(a,A)
     CheckPHRepresentation:the vector and the matrix have different sizes!
     0
    >>> a=[0.2 0.7];
    >>> A=[-1 1; 1 -2];
    >>> CheckPHRepresentation(a,A)
     1
     
     For Python/Numpy:
     
     >>> a=ml.matrix([[0.2]])
     >>> A=ml.matrix([[-1, 1],[1,-2]])
     >>> print(CheckPHRepresentation(a,A))
     CheckPHRepresentation: The vector and the matrix have different sizes!
     False
     >>> a=ml.matrix([[0.2, 0.7]])
     >>> A=ml.matrix([[-1, 1],[1, -2]])
     >>> print(CheckPHRepresentation(a,A))
     True

