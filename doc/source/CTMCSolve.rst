butools.mc.CTMCSolve
====================

.. currentmodule:: butools.mc

.. np:function:: CTMCSolve

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pi = CTMCSolve(Q, prec)`
        * - Mathematica:
          - :code:`pi = CTMCSolve[Q, prec]`
        * - Python/Numpy:
          - :code:`pi = CTMCSolve(Q, prec)`
    
    Computes the stationary solution of a continuous time 
    Markov chain.
    
    Parameters
    ----------
    Q : matrix, shape (M,M)
        The generator matrix of the Markov chain
    prec : double, optional
        Numerical precision for checking whether Q is a 
        valid generator. The default value is 1e-14.
        
    Returns
    -------
    pi : row vector, shape (1,M)
        The vector that satisfies :math:`\pi\, Q = 0, \sum_i \pi_i=1`

    Notes
    -----
    The procedure raises an exception if :code:`butools.checkInput` 
    is set to :code:`true` and :func:`CheckGenerator(Q)` fails.

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
    
    >>> q={{-0.9, 0.5, 0.4}, {0.9, -0.9, 0}, {0.3, 0.3, -0.6}}
    >>> pi=CTMCSolve[q]
    
    For Python/Numpy:

    >>> Q=[[-0.9, 0.5, 0.4], [0.9, -0.9, 0], [0.3, 0.3, -0.6]]
    >>> CTMCSolve(Q)
    matrix([[ 0.40909091,  0.31818182,  0.27272727]])
    
