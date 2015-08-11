butools.mc.DRPSolve
===================

.. currentmodule:: butools.mc

.. np:function:: DRPSolve

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pi = DRPSolve(Q)`
        * - Mathematica:
          - :code:`pi = DRPSolve[Q]`
        * - Python/Numpy:
          - :code:`pi = DRPSolve(Q)`
    
    Computes the stationary solution of a discrete time 
    Markov chain.
    
    Parameters
    ----------
    P : matrix, shape (M,M)
        The matrix parameter of the rational process
        
    Returns
    -------
    pi : row vector, shape (1,M)
        The vector that satisfies 
        :math:`\pi\, P = \pi, \sum_i \pi_i=1`

    Notes
    -----
    Discrete time rational processes are like discrete time 
    Markov chains, but the P matrix does not have to pass 
    the :func:`CheckProbMatrix` test (but the rowsums still 
    have to be ones).

    Examples
    ========
    For Matlab:

    >>> Q = [-0.9, 0.5, 1.4; 0.9, -0.9, 1; 0.3, 1.3, -0.6];
    >>> ret = DRPSolve(Q);
    >>> disp(ret);
          0.23138       0.3484      0.42021
    >>> disp(ret*Q - ret);
       1.6653e-16  -5.5511e-17  -1.1102e-16

    For Mathematica:

    >>> Q = {{-0.9, 0.5, 1.4},{0.9, -0.9, 1},{0.3, 1.3, -0.6}};
    >>> ret = DRPSolve[Q];
    >>> Print[ret];
    {0.2313829787234043, 0.34840425531914887, 0.4202127659574468}
    >>> Print[ret.Q - ret];
    {-1.3877787807814457*^-16, 5.551115123125783*^-17, 5.551115123125783*^-17}

    For Python/Numpy:

    >>> Q = ml.matrix([[-0.9, 0.5, 1.4],[0.9, -0.9, 1],[0.3, 1.3, -0.6]])
    >>> ret = DRPSolve(Q)
    >>> print(ret)
    [[ 0.23138  0.3484   0.42021]]
    >>> print(ret*Q - ret)
    [[ -1.38778e-16   5.55112e-17   5.55112e-17]]

