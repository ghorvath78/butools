butools.queues.MAPMAP1QLD
=========================

.. currentmodule:: butools.queues

.. np:function:: MAPMAP1QLD

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = MAPMAP1QLD(D0, D1, S0, S1, transToPH, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = MAPMAP1QLD[D0, D1, S0, S1, transToPH, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = MAPMAP1QLD(D0, D1, S0, S1, transToPH, prec)`

    Returns the matrix-geometric distribution of the queue
    length of a MAP/MAP/1 queue.

    In a MAP/MAP/1 queue both the arrival and the service
    processes are characterized by Markovian arrival 
    processes.

    Parameters
    ----------
    D0 : matrix, shape(N,N)
        The transitions of the arrival MAP not accompanied by
        job arrivals
    D1 : matrix, shape(N,N)
        The transitions of the arrival MAP accompanied by
        job arrivals
    S0 : matrix, shape(N,N)
        The transitions of the service MAP not accompanied by
        job service
    S1 : matrix, shape(N,N)
        The transitions of the service MAP accompanied by
        job service
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
    >>> S0 = [-10, 4; 0, -7];
    >>> S1 = [5, 1; 4, 3];
    >>> [alpha,A] = MAPMAP1QLD(D0,D1,S0,S1);
    >>> MomentsFromMG(alpha,A,5)
          0.54864        1.306        4.357       19.193       105.39
    >>> x=(0:10)';
    >>> y=PmfFromMG(alpha,A,x);
    >>> plot(x,y);
    >>> [alpha,A] = MAPMAP1QLD(D0,D1,S0,S1,true);
    >>> CheckDPHRepresentation(alpha,A)
         1    
    
