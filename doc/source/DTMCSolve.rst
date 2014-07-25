butools.mc.DTMCSolve
====================

.. currentmodule:: butools.mc

.. np:function:: DTMCSolve

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pi = DTMCSolve(Q, prec)`
        * - Mathematica:
          - :code:`pi = DTMCSolve[Q, prec]`
        * - Python/Numpy:
          - :code:`pi = DTMCSolve(Q, prec)`
    
    Computes the stationary solution of a discrete time 
    Markov chain.
    
    Parameters
    ----------
    P : matrix, shape (M,M)
        The transition probability matrix of the Markov 
        chain
    prec : double, optional
        Numerical precision for checking whether P is a 
        valid generator. The default value is 1e-14.
        
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
    
    >>> q={{0.1, 0.5, 0.4}, {0.9, 0.1, 0}, {0.3, 0.3, 0.4}}
    >>> pi=DTMCSolve[q]
    
    For Python/Numpy:

    >>> q=[[0.1, 0.5, 0.4], [0.9, 0.1, 0], [0.3, 0.3, 0.4]]
    >>> pi=DTMCSolve(q)
        
