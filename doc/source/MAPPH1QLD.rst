butools.queues.MAPPH1QLD
========================

.. currentmodule:: butools.queues

.. np:function:: MAPPH1QLD

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = MAPPH1QLD(D0, D1, sigma, S, transToPH, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = MAPPH1QLD[D0, D1, sigma, S, transToPH, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = MAPPH1QLD(D0, D1, sigma, S, transToPH, prec)`

    Returns the matrix-geometric distribution of the queue
    length of a MAP/PH/1 queue.

    In a MAP/PH/1 queue the arrival processes is given by a 
    Markovian arrival processes, and the service time follows
    a phase-type distribution.

    Parameters
    ----------
    D0 : matrix, shape(N,N)
        The transitions of the arrival MAP not accompanied by
        job arrivals
    D1 : matrix, shape(N,N)
        The transitions of the arrival MAP accompanied by
        job arrivals
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

    >>> D0 = [-8, 1, 2; 0, -6, 4; 3 0 -3];
    >>> D1 = [4, 1, 0; 0, 2, 0; 0, 0, 0];
    >>> sigma = [0.2, 0.7, 0.1];
    >>> S = [-10, 4, 0; 5, -7, 2; 1, 2, -8];
    >>> [alpha,A] = MAPPH1QLD(D0,D1,sigma,S);
    >>> MomentsFromMG(alpha,A,5)
           3.8041       33.872       450.13       7974.7    1.766e+05
    >>> x=(0:15)';
    >>> y=PmfFromMG(alpha,A,x);
    >>> plot(x,y);
    >>> [alpha,A] = MAPPH1QLD(D0,D1,sigma,S,true);
    >>> CheckDPHRepresentation(alpha,A)
         1    
    
