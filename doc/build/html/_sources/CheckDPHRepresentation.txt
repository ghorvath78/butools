butools.dph.CheckDPHRepresentation
==================================

.. currentmodule:: butools.dph

.. np:function:: CheckDPHRepresentation

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`r = CheckDPHRepresentation(alpha, A, prec)`
        * - Mathematica:
          - :code:`r = CheckDPHRepresentation[alpha, A, prec]`
        * - Python/Numpy:
          - :code:`r = CheckDPHRepresentation(alpha, A, prec)`

    Checks if the given vector and matrix define a valid 
    discrete phase-type representation.

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
        A is substochastic, and they have the same size.

    Examples
    ========
    For Matlab:

    >>> a = [0.48, 0.08, 0.26, 0.18];
    >>> A = [0, 0.08, 0.08, 0.8; 0.55, 0, 0.24, 0.19; 0.06, 0.03, 0, 0.001; 0.23, 0.005, 0.2, 0.53];
    >>> flag = CheckDPHRepresentation(a,A);
    >>> disp(flag);
         1
    >>> a = [0.48, 0.08];
    >>> A = [0, 0.08; 0.55, 0.5];
    >>> flag = CheckDPHRepresentation(a,A);
    CheckProbMatrix: The rowsum of the matrix (transient) is not less or equal than 1 (precision: 1e-12)!
    >>> disp(flag);
         0

    For Mathematica:

    >>> a = {0.48, 0.08, 0.26, 0.18};
    >>> A = {{0, 0.08, 0.08, 0.8},{0.55, 0, 0.24, 0.19},{0.06, 0.03, 0, 0.001},{0.23, 0.005, 0.2, 0.53}};
    >>> flag = CheckDPHRepresentation[a,A];
    >>> Print[flag];
    True
    >>> a = {0.48, 0.08};
    >>> A = {{0, 0.08},{0.55, 0.5}};
    >>> flag = CheckDPHRepresentation[a,A];
    "CheckProbMatrix: A rowsum of the transient matrix is not less or equal than 1!"
    >>> Print[flag];
    False

    For Python/Numpy:

    >>> a = ml.matrix([[0.48, 0.08, 0.26, 0.18]])
    >>> A = ml.matrix([[0, 0.08, 0.08, 0.8],[0.55, 0, 0.24, 0.19],[0.06, 0.03, 0, 0.001],[0.23, 0.005, 0.2, 0.53]])
    >>> flag = CheckDPHRepresentation(a,A)
    >>> print(flag)
    True
    >>> a = ml.matrix([[0.48, 0.08]])
    >>> A = ml.matrix([[0, 0.08],[0.55, 0.5]])
    >>> flag = CheckDPHRepresentation(a,A)
    CheckProbMatrix: The rowsum of the matrix (transient) is not less or equal than 1 (precision: {0})! 1e-12
    >>> print(flag)
    False

