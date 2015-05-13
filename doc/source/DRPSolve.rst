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
    --------
    For Matlab:
    
    >>> Q=[-0.9 0.5 1.4; 0.9 -0.9 1; 0.3 1.3 -0.6];
    >>> pi=DRPSolve(Q)
    [0.23138 0.3484 0.42021]
    >>> sum(pi)
    1
    >>> pi*Q - pi
    [1.6653e-16 -5.5511e-17 -1.1102e-16]    
    
    For Mathematica:
    
    >>> q={{-0.9, 0.5, 1.4}, {0.9, -0.9, 1}, {0.3, 1.3, -0.6}};
    >>> pi=DRPSolve[q]
    {0.231383, 0.348404, 0.420213}
    
    For Python/Numpy:

    >>> Q=[[-0.9, 0.5, 1.4], [0.9, -0.9, 1], [0.3, 1.3, -0.6]]
    >>> DRPSolve(Q)
    matrix([[ 0.23138298,  0.34840426,  0.42021277]])

