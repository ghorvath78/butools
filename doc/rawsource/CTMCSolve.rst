butools.mc.CTMCSolve
====================

.. currentmodule:: butools.mc

.. np:function:: CTMCSolve

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pi = CTMCSolve(Q)`
        * - Mathematica:
          - :code:`pi = CTMCSolve[Q]`
        * - Python/Numpy:
          - :code:`pi = CTMCSolve(Q)`
    
    Computes the stationary solution of a continuous time 
    Markov chain.
    
    Parameters
    ----------
    Q : matrix, shape (M,M)
        The generator matrix of the Markov chain
        
    Returns
    -------
    pi : row vector, shape (1,M)
        The vector that satisfies :math:`\pi\, Q = 0, \sum_i \pi_i=1`

    Notes
    -----
    The procedure raises an exception if :code:`checkInput` 
    is set to :code:`true` and :func:`CheckGenerator` (Q) fails.

    Examples
    --------
    For Matlab:
    
    >>> Q=[-0.9 0.5 0.4; 0.9 -0.9 0; 0.3 0.3 -0.6]
    >>> pi=CTMCSolve(Q)
    [0.4091 0.3182 0.2727]
    >>> sum(pi)
    1
    >>> pi*Q
    [-1.1102e-16 1.3878e-17 8.3267e-17]
    
    For Mathematica:
    
    >>> q={{-0.9, 0.5, 0.4}, {0.9, -0.9, 0}, {0.3, 0.3, -0.6}};
    >>> pi=CTMCSolve[q]
    {0.409091, 0.318182, 0.272727}
    >>> q={{-9/10, 5/10, 4/10}, {9/10, -9/10, 0}, {3/10, 3/10, -6/10}};
    >>> pi=CTMCSolve[q]
    {9/22, 7/22, 3/11}
    
    For Python/Numpy:

    >>> Q=[[-0.9, 0.5, 0.4], [0.9, -0.9, 0], [0.3, 0.3, -0.6]]
    >>> CTMCSolve(Q)
    matrix([[ 0.40909091,  0.31818182,  0.27272727]])
    
