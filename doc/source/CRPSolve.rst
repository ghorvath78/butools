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
    --------
    For Matlab:
    
    >>> Q=[-4.3 3.5 0.8; -8.4 6.5 1.9; 17.3 -12.7 -4.6];
    >>> pi=CRPSolve(Q)
    [-3.5617 3.6667 0.89506]
    >>> sum(pi)
    1
    >>> pi*Q
    [3.5527e-15 -1.7764e-15 0]   
    
    For Mathematica:
    
    >>> q={{-4.3, 3.5, 0.8}, {-8.4, 6.5, 1.9}, {17.3, -12.7, -4.6}};
    >>> pi=CRPSolve[q]
    {-3.56173, 3.66667, 0.895062}
    
    For Python/Numpy:

    >>> Q=[[-4.3, 3.5, 0.8],[-8.4, 6.5, 1.9],[17.3, -12.7, -4.6]]
    >>> CRPSolve(Q)
    matrix([[-3.5617284 ,  3.66666667,  0.89506173]])

