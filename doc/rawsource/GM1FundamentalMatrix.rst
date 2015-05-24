butools.mam.GM1FundamentalMatrix
================================

.. currentmodule:: butools.mam

.. np:function:: GM1FundamentalMatrix

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`R = GM1FundamentalMatrix(A, precision, maxNumIt, method)`
        * - Mathematica:
          - :code:`R = GM1FundamentalMatrix[A, precision, maxNumIt, method]`
        * - Python/Numpy:
          - :code:`R = GM1FundamentalMatrix(A, precision, maxNumIt, method)`

    Returns matrix R corresponding to the G/M/1 type Markov
    chain given by matrices A.
    
    Matrix R is the minimal non-negative solution of the 
    following matrix equation:
    
    .. math::
        R = A_0 + R A_1 + R^2 A_2 + R^3 A_3 + \dots.
    
    The implementation is based on [1]_, please cite it if
    you use this method.
    
    Parameters
    ----------
    A : matrix, shape (N,M*N)
        Matrix blocks of the G/M/1 type generator from
        0 to M-1, concatenated horizontally.
    precision : double, optional
        Matrix R is computed iteratively up to this
        precision. The default value is 1e-14
    maxNumIt : int, optional
        The maximal number of iterations. The default value
        is 50.
    method : {"CR", "RR", "NI", "FI", "IS"}, optional
        The method used to solve the matrix-quadratic
        equation (CR: cyclic reduction, RR: Ramaswami
        reduction, NI: Newton iteration, FI: functional
        iteration, IS: invariant subspace method). The 
        default is "CR".
    
    Returns
    -------
    R : matrix, shape (N,N)
        The R matrix of the G/M/1 type Markov chain.
    
    References
    ----------
    .. [1] Bini, D. A., Meini, B., Steffé, S., Van Houdt,
           B. (2006, October). Structured Markov chains 
           solver: software tools. In Proceeding from the
           2006 workshop on Tools for solving structured 
           Markov chains (p. 14). ACM.

    Examples
    --------
    For Matlab:
    
    >>> A0 = [0.1, 0; 0, 0.1];
    >>> A1 = [0, 0.2; 0, 0.2];
    >>> A2 = [0, 0.1; 0, 0];
    >>> A3 = [0.3, 0.2; 0.3, 0.2];
    >>> A4 = [0, 0.1; 0.2, 0];
    >>> A = [A0,A1,A2,A3,A4];
    >>> R = GM1FundamentalMatrix(A);
          0.10065     0.026961
       0.00065531      0.12569

    For Python/Numpy:
       
    >>> A0 = ml.matrix([[0.1, 0],[0, 0.1]])
    >>> A1 = ml.matrix([[0, 0.2],[0, 0.2]])
    >>> A2 = ml.matrix([[0, 0.1],[0, 0]])
    >>> A3 = ml.matrix([[0.3, 0.2],[0.3, 0.2]])
    >>> A4 = ml.matrix([[0, 0.1],[0.2, 0]])
    >>> A = (A0,A1,A2,A3,A4)
    >>> R = GM1FundamentalMatrix (A)
    >>> print(R)
    [[ 0.1006515   0.02696092]
     [ 0.00065531  0.1256871 ]]

