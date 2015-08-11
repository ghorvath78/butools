butools.fitting.MAPFromTrace
============================

.. currentmodule:: butools.fitting

.. np:function:: MAPFromTrace

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[D0, D1, logli] = MAPFromTrace(trace, orders, maxIter, stopCond, initial, result)`
        * - Mathematica:
          - :code:`{D0, D1, logli} = MAPFromTrace[trace, orders, maxIter, stopCond, initial, result]`
        * - Python/Numpy:
          - :code:`D0, D1, logli = MAPFromTrace(trace, orders, maxIter, stopCond, initial, result)`

    Performs MAP fitting using the EM algorithm (ErCHMM, 
    [1]_, [2]_).
    
    Parameters
    ----------
    trace : column vector, length K
        The samples of the trace
    orders : list of int, length(N), or int
        The length of the list determines the number of 
        Erlang branches to use in the fitting method.
        The entries of the list are the orders of the 
        Erlang distributions. If this parameter is a 
        single integer, all possible branch number - order
        combinations are tested where the total number of 
        states is "orders".
    maxIter : int, optional
        Maximum number of iterations. The default value is
        200
    stopCond : double, optional
        The algorithm stops if the relative improvement of
        the log likelihood falls below stopCond. The 
        default value is 1e-7
    initial : tuple of a vector and a matrix, shape(N,N), optional
        The rate parameters of the Erlang distributions 
        and the branch transition probability matrix to be
        used initially. If not given, a default initial 
        guess is determined and the algorithm starts from 
        there.
    result : {"vecmat", "matmat"}, optional
        The result can be returned two ways. If "matmat" is
        selected, the result is returned in the classical
        representation of MAPs, thus the D0 and D1 matrices.
        If "vecmat" is selected, the rate parameters of the
        Erlang branches and the branch transition probability
        matrix are returned. The default value is "matmat"

    Returns
    -------
    (D0, D1) : tuple of matrix, shape (M,M) and matrix, shape (M,M)
        If the "matmat" result format is chosen, the function
        returns the D0 and D1 matrices of the MAP
    (lambda, P) : tuple of vector, length N and matrix, shape (M,M)
        If the "vecmat" result format is chosen, the function
        returns the vector of the Erlang rate parameters of 
        the branches and the branch transition probability 
        matrix
    logli : double
        The log-likelihood divided by the trace length
        
    Notes
    -----
    This procedure is quite slow in the supported 
    mathematical frameworks. If the maximum speed is
    needed, please use the multi-core optimized c++
    implementation called SPEM-FIT_.

    .. _SPEM-FIT: https://bitbucket.org/ghorvath78/spemfit

    References
    ----------
    .. [1] Okamura, Hiroyuki, and Tadashi Dohi. Faster 
           maximum likelihood estimation algorithms for 
           Markovian arrival processes. Quantitative 
           Evaluation of Systems, 2009. QEST'09. Sixth 
           International Conference on the. IEEE, 2009.
    
    .. [2] Horváth, Gábor, and Hiroyuki Okamura. A Fast EM
           Algorithm for Fitting Marked Markovian Arrival 
           Processes with a New Special Structure. Computer
           Performance Engineering. Springer Berlin 
           Heidelberg, 2013. 119-133.
    
    Examples
    --------    
    For Matlab:
    
    >>> tr = dlmread('trace.txt');
    >>> [D0,D1] = MAPFromTrace(tr(1:10000),5);
    >>> D0
          -83.429            0            0            0            0
                0      -718.68            0            0            0
                0            0      -1026.2       1026.2            0
                0            0            0      -1026.2       1026.2
                0            0            0            0      -1026.2
    >>> D1
           54.149       4.9019       24.379            0            0
           3.3915       665.85       49.439            0            0
                0            0            0            0            0
                0            0            0            0            0
           42.647       96.944       886.57            0            0
    
    For Python/Numpy:
    
    >>> tr = np.loadtxt('trace.txt')
    >>> D0,D1 = MAPFromTrace(tr[0:10000],5)
    >>> print(D0)
    [[  -83.42943389     0.             0.             0.             0.        ]
     [    0.          -718.67798921     0.             0.             0.        ]
     [    0.             0.         -1026.16062945  1026.16062945     0.        ]
     [    0.             0.             0.         -1026.16062945  1026.16062945]
     [    0.             0.             0.             0.         -1026.16062945]]
    >>> print(D1)
    [[  54.14857309    4.90186422   24.37899658    0.            0.        ]
     [   3.39151863  665.84734646   49.43912412    0.            0.        ]
     [   0.            0.            0.            0.            0.        ]
     [   0.            0.            0.            0.            0.        ]
     [  42.64730202   96.94396094  886.5693665     0.            0.        ]]

