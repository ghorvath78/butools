butools.queues.PHPH1QLD
=======================

.. currentmodule:: butools.queues

.. np:function:: PHPH1QLD

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = PHPH1QLD(delta, D, sigma, S, transToPH, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = PHPH1QLD[delta, D, sigma, S, transToPH, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = PHPH1QLD(delta, D, sigma, S, transToPH, prec)`

    Returns the matrix-geometric distribution of the queue
    length of a PH/PH/1 queue.

    In a PH/PH/1 queue both the inter-arrival times and the 
    service times are given by phase-type distributions.

    Parameters
    ----------
    delta : matrix, shape(1,N)
        The initial probability vector of the phase-type
        distribution generating the inter-arrival times
    D : matrix, shape(N,N)
        The transient generator of the phase-type 
        distribution generating the inter-arrival times
    sigma : matrix, shape(1,N)
        The initial probability vector of the phase-type
        distribution corresponging to the service time
    S : matrix, shape(N,N)
        The transient generator of the phase-type 
        distribution corresponging to the service time
    transToPH : bool, optional
        If true, the result is transformed to a discrete 
        phase-type representation. The default value is false
    prec : double, optional
        Numerical precision to check if the input is valid and
        it is also used as a stopping condition when solving
        the matrix-quadratic equation

    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial vector of the matrix-geometric 
        distribution of the number of jobs in the system
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-geometric
        distribution of the number of jobs in the system

    Notes
    -----
    While the transformation of the results to a discrete 
    phase-type representation is always possible 
    theoretically, it may introduce numerical problems in 
    some cases.

    Examples
    ========
    For MATLAB:

    >>> delta = [0.5, 0.1, 0.4];
    >>> D = [-8, 1, 2; 0, -6, 4; 3 0 -3];
    >>> sigma = [0.2, 0.7, 0.1];
    >>> S = [-10, 4, 0; 5, -7, 2; 1, 2, -8];
    >>> [alpha,A] = PHPH1QLD(delta,D,sigma,S);
    >>> MomentsFromMG(alpha,A,5)
           2.0439       10.554       80.619       820.69        10443
    >>> x=(0:10)';
    >>> y=PmfFromMG(alpha,A,x);
    >>> plot(x,y);
    >>> [alpha,A] = PHPH1QLD(delta,D,sigma,S,true);
    >>> CheckDPHRepresentation(alpha,A)
         1    
    
