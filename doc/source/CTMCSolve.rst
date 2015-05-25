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
    The procedure raises an exception if :code:`butools.checkInput` 
    is set to :code:`true` and :func:`CheckGenerator(Q)` fails.

    Examples
    ========
    For Matlab:

    >>> Q = [-0.9, 0.5, 0.4; 0.9, -0.9, 0; 0.3, 0.3, -0.6];
    >>> ret = CTMCSolve(Q);
    >>> disp(ret);
          0.40909      0.31818      0.27273
    >>> disp(ret*Q);
      -1.1102e-16   1.3878e-17   8.3267e-17

    For Mathematica:

    >>> Q = {{-0.9, 0.5, 0.4},{0.9, -0.9, 0},{0.3, 0.3, -0.6}};
    >>> ret = CTMCSolve[Q];
    >>> Print[ret];
    {0.4090909090909091, 0.3181818181818182, 0.2727272727272727}
    >>> Print[ret.Q];
    {-4.163336342344337*^-17, -1.3877787807814457*^-17, 5.551115123125783*^-17}

    For Python/Numpy:

    >>> Q = ml.matrix([[-0.9, 0.5, 0.4],[0.9, -0.9, 0],[0.3, 0.3, -0.6]])
    >>> ret = CTMCSolve(Q)
    >>> print(ret)
    [[ 0.40909  0.31818  0.27273]]
    >>> print(ret*Q)
    [[ -4.16334e-17  -1.38778e-17   5.55112e-17]]

