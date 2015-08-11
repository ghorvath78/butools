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
    ========
    For Matlab:

    >>> Q = [0.1, 0.5, 0.4; 0.9, 0.1, 0; 0.3, 0.3, 0.4];
    >>> ret = DTMCSolve(Q);
    >>> disp(ret);
          0.40909      0.31818      0.27273
    >>> disp(ret*Q - ret);
      -5.5511e-17   5.5511e-17   5.5511e-17

    For Mathematica:

    >>> Q = {{0.1, 0.5, 0.4},{0.9, 0.1, 0},{0.3, 0.3, 0.4}};
    >>> ret = DTMCSolve[Q];
    >>> Print[ret];
    {0.4090909090909091, 0.3181818181818182, 0.2727272727272727}
    >>> Print[ret.Q - ret];
    {-5.551115123125783*^-17, 0., 5.551115123125783*^-17}

    For Python/Numpy:

    >>> Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, 0.3, 0.4]])
    >>> ret = DTMCSolve(Q)
    >>> print(ret)
    [[ 0.40909  0.31818  0.27273]]
    >>> print(ret*Q - ret)
    [[ -5.55112e-17   0.00000e+00   5.55112e-17]]

