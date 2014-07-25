butools.fitting.PHFromTrace
===========================

.. currentmodule:: butools.fitting

.. np:function:: PHFromTrace

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A, logli] = PHFromTrace(trace, orders, maxIter, stopCond, initial, result)`
        * - Mathematica:
          - :code:`{alpha, A, logli} = PHFromTrace[trace, orders, maxIter, stopCond, initial, result]`
        * - Python/Numpy:
          - :code:`alpha, A, logli = PHFromTrace(trace, orders, maxIter, stopCond, initial, result)`

    Performs PH distribution fitting using the EM algorithm
    (G-FIT, [1]_).
    
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
    initial : tuple of two vectors, optional
        The initial values of the branch probabilities and
        rate parameters is given by this tuple. If not 
        given, a default initial guess is determined and 
        the algorithm starts from there.
    result : {"vecmat", "vecvec"}, optional
        The result can be returned two ways. If "vecmat" is
        selected, the result is returned in the classical
        representation of phase-type distributions, thus the
        initial vector and the generator matrix. 
        If "vecvec" is selected, two vectors are returned, 
        one holds the branch probabilities, and the second
        holds the rate parameters of the Erlang branches.
        The default value is "vecmat"

    Returns
    -------
    (alpha, A) : tuple of matrix, shape (1,M) and matrix, shape (M,M)
        If the "vecmat" result format is chosen, the function
        returns the initial probability vector and the
        generator matrix of the phase type distribution.
    (pi, lambda) : tuple of vector, length N and vector, length N
        If the "vecvec" result format is chosen, the function
        returns the vector of branch probabilities and the
        vector of branch rates in a tuple.
    logli : double
        The log-likelihood divided by the trace length
        
    Notes
    -----
    This procedure is quite fast in the supported 
    mathematical frameworks. If the maximum speed is
    needed, please use the multi-core optimized c++
    implementation called SPEM-FIT_.

    .. _SPEM-FIT: https://bitbucket.org/ghorvath78/spemfit

    References
    ----------
    .. [1] Thummler, Axel, Peter Buchholz, and Miklós Telek.
           A novel approach for fitting probability 
           distributions to real trace data with the EM 
           algorithm. Dependable Systems and Networks, 2005.

    Examples
    --------    
    For Matlab:
    
    >>> tr = dlmread('trace.txt');
    >>> [alpha,A] = PHFromTrace(tr,5);
    >>> alpha
         0.065027      0.85788            0     0.077088            0
    >>> A
          -63.308            0            0            0            0
                0      -815.72       815.72            0            0
                0            0      -815.72            0            0
                0            0            0       -12563        12563
                0            0            0            0       -12563
    
    
