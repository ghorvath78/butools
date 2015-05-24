butools.ph.CheckMERepresentation
================================

.. currentmodule:: butools.ph

.. np:function:: CheckMERepresentation

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckMERepresentation(alpha, A, prec)`
        * - Mathematica:
          - :code:`r = CheckMERepresentation[alpha, A, prec]`
        * - Python/Numpy:
          - :code:`r = CheckMERepresentation(alpha, A, prec)`

    Checks if the given vector and matrix define a valid matrix-
    exponential representation.

    Parameters
    ----------
    alpha : matrix, shape (1,M)
        Initial vector of the matrix-exponential distribution 
        to check
    A : matrix, shape (M,M)
        Matrix parameter of the matrix-exponential distribution
        to check
    prec : double, optional
        Numerical precision. The default value is 1e-14.

    Returns
    -------
    r : bool
        True, if the matrix is a square matrix, the vector and 
        the matrix have the same size, the dominant eigenvalue
        is negative and real
    
    Notes
    -----
    This procedure does not check the positivity of the density!
    Call 'CheckMEPositiveDensity' if it is needed, but keep in
    mind that it can be time-consuming, while this procedure
    is fast.

    Examples
    --------
    For Matlab:
    
    >>> a=[-0.2 0.2];
    >>> A=[1 -1; 1 -2];
    >>> CheckMERepresentation(a,A)
     CheckMERepresentation: There is an eigenvalue of the matrix with non-negative real part (at precision 1e-14)!
     0   
    >>> a=[-0.2 0.4 0.8];
    >>> A=[-2 0 3; 0 -1 1; 0 -1 -1];
    >>> CheckMERepresentation(a,A)
     CheckMERepresentation: The dominant eigenvalue of the matrix is not real!
     0      
    >>> a = [0.2, 0.3, 0.5];
    >>> A = [-1,0,0;0,-3,2;0,-2,-3];
    >>> CheckMERepresentation(a,A)
     1

    For Python/Numpy:
    
    >>> a=ml.matrix([[-0.2, 0.2]])
    >>> A=ml.matrix([[1, -1],[1, -2]])
    >>> print(CheckMERepresentation(a,A))
    CheckMERepresentation: There is an eigenvalue of the matrix with non-negative real part!
    False
    >>> a=ml.matrix([[-0.2, 0.4, 0.8]])
    >>> A=ml.matrix([[-2, 0, 3],[0, -1, 1],[0, -1, -1]])
    >>> print(CheckMERepresentation(a,A))
    CheckMERepresentation: The dominant eigenvalue of the matrix is not real!
    False
    >>> a = ml.matrix([[0.2, 0.3, 0.5]])
    >>> A = ml.matrix([[-1,0,0],[0,-3,2],[0,-2,-3]])
    >>> print(CheckMERepresentation(a,A))
    True

