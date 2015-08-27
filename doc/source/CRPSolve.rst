butools.mc.CRPSolve
===================

.. currentmodule:: butools.mc

.. np:function:: CRPSolve

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pi = CRPSolve(Q)`
        * - Mathematica:
          - :code:`pi = CRPSolve[Q]`
        * - Python/Numpy:
          - :code:`pi = CRPSolve(Q)`
    
    Computes the stationary solution of a continuous time 
    rational process (CRP).
    
    Parameters
    ----------
    Q : matrix, shape (M,M)
        The generator matrix of the rational process
        
    Returns
    -------
    pi : row vector, shape (1,M)
        The vector that satisfies 
        :math:`\pi\, Q = 0, \sum_i \pi_i=1`

    Notes
    -----
    Continuous time rational processes are like continuous 
    time Markov chains, but the generator does not have to 
    pass the :func:`CheckGenerator` test (but the rowsums 
    still have to be zeros).

    Examples
    ========
    For Matlab:

    >>> Q = [-4.3, 3.5, 0.8; -8.4, 6.5, 1.9; 17.3, -12.7, -4.6];
    >>> ret = CRPSolve(Q);
    >>> disp(ret);
          -3.5617       3.6667      0.89506
    >>> disp(ret*Q);
       3.5527e-15  -1.7764e-15            0

    For Mathematica:

    >>> Q = {{-4.3, 3.5, 0.8},{-8.4, 6.5, 1.9},{17.3, -12.7, -4.6}};
    >>> ret = CRPSolve[Q];
    >>> Print[ret];
    {-3.5617283950617336, 3.66666666666667, 0.8950617283950623}
    >>> Print[ret.Q];
    {1.7763568394002505*^-15, 0., 0.}

    For Python/Numpy:

    >>> Q = ml.matrix([[-4.3, 3.5, 0.8],[-8.4, 6.5, 1.9],[17.3, -12.7, -4.6]])
    >>> ret = CRPSolve(Q)
    >>> print(ret)
    [[-3.56173  3.66667  0.89506]]
    >>> print(ret*Q)
    [[  1.77636e-15   0.00000e+00   0.00000e+00]]

