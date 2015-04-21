butools.mam.QBDFundamentalMatrices
==================================

.. currentmodule:: butools.mam

.. np:function:: QBDFundamentalMatrices

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`M = QBDFundamentalMatrices(B, L, F, matrices, precision, maxNumIt, method)`
        * - Mathematica:
          - :code:`M = QBDFundamentalMatrices[xxx]`
        * - Python/Numpy:
          - :code:`M = QBDFundamentalMatrices(B, L, F, matrices, precision, maxNumIt, method)`

    Returns the fundamental matrices corresponding to the
    given QBD. Matrices R, G and U can be returned
    depending on the "matrices" parameter.
    
    The implementation is based on [1]_, please cite it if
    you use this method.
    
    Parameters
    ----------
    B : matrix, shape (N,N)
        The matrix corresponding to backward transitions
    L : matrix, shape (N,N)
        The matrix corresponding to local transitions
    F : matrix, shape (N,N)
        The matrix corresponding to forward transitions
    matrices : string
        Specifies which matrices are required. 'R' means 
        that only matrix 'R' is needed. 'UG' means that
        matrices U and G are needed. Any combinations of
        'R', 'G' and 'U' are allowed, in any order.
    precision : double, optional
        The matrices are computed iteratively up to this
        precision. The default value is 1e-14
    maxNumIt : int, optional
        The maximal number of iterations. The default value
        is 50.
    method : {"CR", "LR", "NI", "FI", "IS"}, optional
        The method used to solve the matrix-quadratic
        equation (CR: cyclic reduction, LR: logarithmic
        reduction, NI: Newton iteration, FI: functional
        iteration, IS: invariant subspace method). The 
        default is "CR".
    
    Returns
    -------
    M : list of matrices
        The list of calculated matrices in the order as
        requested in the 'matrices' parameter.
    
    Notes
    -----
    Discrete and continuous QBDs are both supported, the
    procedure auto-detects it based on the diagonal entries
    of matrix L.

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
    
    >>> B = [0,0;3,4];
    >>> L = [-6,5;3,-12];
    >>> F = [1,0;2,0];
    >>> [R, G, U] = QBDFundamentalMatrices (B,L,F,'RGU');
    >>> R
          0.27839      0.14286
          0.55678      0.28571
    >>> G
          0.42857      0.57143
          0.42857      0.57143
    >>> U
          -5.5714       5.5714
           3.8571      -10.857

    For Python/Numpy:
    
    >>> B = ml.matrix([[0,0],[3,4]])
    >>> L = ml.matrix([[-6,5],[3,-12]])
    >>> F = ml.matrix([[1,0],[2,0]])
    >>> R, G, U = QBDFundamentalMatrices (B,L,F,"RGU")
    >>> print(R)
    [[ 0.27838828  0.14285714]
     [ 0.55677656  0.28571429]]
    >>> print(G)
    [[ 0.42857143  0.57142857]
     [ 0.42857143  0.57142857]]
    >>> print(U)
    [[ -5.57142857   5.57142857]
     [  3.85714286 -10.85714286]]

