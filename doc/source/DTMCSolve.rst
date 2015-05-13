butools.mc.DTMCSolve
====================

.. currentmodule:: butools.mc

.. np:function:: DTMCSolve

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pi = DTMCSolve(Q)`
        * - Mathematica:
          - :code:`pi = DTMCSolve[Q]`
        * - Python/Numpy:
          - :code:`pi = DTMCSolve(Q)`
    
    Computes the stationary solution of a discrete time 
    Markov chain.
    
    Parameters
    ----------
    P : matrix, shape (M,M)
        The transition probability matrix of the Markov 
        chain
        
    Returns
    -------
    pi : row vector, shape (1,M)
        The vector that satisfies :math:`\pi\, P = \pi, \sum_i \pi_i=1`

    Notes
    -----
    The procedure raises an exception if :code:`butools.checkInput` 
    is set to :code:`true` and :func:`CheckProbMatrix(P)` fails.

    Examples
    --------
    For Matlab:
    
    >>> Q=[0.1 0.5 0.4; 0.9 0.1 0; 0.3 0.3 0.4];
    >>> pi=DTMCSolve(Q)
    [0.40909 0.31818 0.27273]
    >>> sum(pi)
    1
    >>> pi*Q - pi
    [-5.5511e-17 5.5511e-17 5.5511e-17]    
    
    For Mathematica:
    
    >>> Q={{0.1, 0.5, 0.4}, {0.9, 0.1, 0}, {0.3, 0.3, 0.4}};
    >>> pi=DTMCSolve[Q]
    {0.409091, 0.318182, 0.272727}
    >>> Q = {{1/10, 5/10, 4/10}, {9/10, 1/10, 0}, {3/10, 3/10, 4/10}};
    >>> pi=DTMCSolve[Q]
    {9/22, 7/22, 3/11}
    
    For Python/Numpy:

    >>> P=[[0.1, 0.5, 0.4], [0.9, 0.1, 0], [0.3, 0.3, 0.4]]
    >>> DTMCSolve(P)
    matrix([[ 0.40909091,  0.31818182,  0.27272727]])
        
