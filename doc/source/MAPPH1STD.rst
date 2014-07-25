butools.queues.MAPPH1STD
========================

.. currentmodule:: butools.queues

.. np:function:: MAPPH1STD

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = MAPPH1STD(D0, D1, sigma, S, transToPH, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = MAPPH1STD[D0, D1, sigma, S, transToPH, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = MAPPH1STD(D0, D1, sigma, S, transToPH, prec)`

    Returns the matrix-exponential distribution of the sojourn
    time of the jobs in a MAP/PH/1 queue.

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
        If true, the result is transformed to a phase-type
        representation. The default value is false
    prec : double, optional
        Numerical precision to check if the input is valid and
        it is also used as a stopping condition when solving
        the matrix-quadratic equation

    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial vector of the matrix-exponential 
        distribution representing the sojourn time of the jobs
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential 
        distribution representing the sojourn time of the jobs
        

    Notes
    -----
    While the transformation of the results to phase-type
    representation is always possible theoretically, it may
    introduce numerical problems in some cases.

    Examples
    ========
    For MATLAB:

    >>> D0 = [-8, 1, 2; 0, -6, 4; 3 0 -3];
    >>> D1 = [4, 1, 0; 0, 2, 0; 0, 0, 0];
    >>> sigma = [0.2, 0.7, 0.1];
    >>> S = [-10, 4, 0; 5, -7, 2; 1, 2, -8];
    >>> [alpha,A] = MAPPH1STD(D0,D1,sigma,S);
    >>> MomentsFromME(alpha,A,5)
           1.7902       6.2919       33.139       232.71       2042.6
    >>> x=(0:0.1:4)';
    >>> y=PdfFromME(alpha,A,x);
    >>> plot(x,y);
    >>> [alpha,A] = MAPPH1STD(D0,D1,sigma,S,true);
    >>> CheckPHRepresentation(alpha,A)
         1    
    
