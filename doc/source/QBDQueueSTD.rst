butools.queues.QBDQueueSTD
==========================

.. currentmodule:: butools.queues

.. np:function:: QBDQueueSTD

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = QBDQueueSTD(B, L, F, L0, transToPH, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = QBDQueueSTD[B, L, F, L0, transToPH, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = QBDQueueSTD(B, L, F, L0, transToPH, prec)`

    Returns the matrix-exponential distribution of the sojourn
    time of the jobs in a QBD queue.

    QBD queues have a background continuous time Markov chain
    with generator Q whose the transitions can be partitioned
    into three sets: transitions accompanied by an arrival
    of a new job (F, forward), transitions accompanied by 
    the service of the current job in the server (B, 
    backward) and internal transitions (L, local). 
    Thus we have Q=B+L+F.

    Parameters
    ----------
    B : matrix, shape(N,N)
        Transitions of the background process accompanied by 
        the service of the current job in the server
    L : matrix, shape(N,N)
        Internal transitions of the background process 
        that do not generate neither arrival nor service
    F : matrix, shape(N,N)
        Transitions of the background process accompanied by 
        an arrival of a new job
    L0 : matrix, shape(N,N)
        Internal transitions of the background process when
        there are no jobs in the queue
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

    >>> B=[6, 1, 0; 0, 4, 1; 2, 0, 0];
    >>> F=[0, 0, 0; 5, 0, 0; 1, 3, 0];
    >>> L=[-12, 3, 2; 0, -14, 4; 3, 1, -10];
    >>> L0=L+B;
    >>> [alpha,A] = QBDQueueSTD(B, L, F, L0);
    >>> MomentsFromME(alpha,A,5)
          0.33414      0.23413      0.24869      0.35297      0.62645
    >>> x=(0:0.05:1)';
    >>> y=PdfFromME(alpha,A,x);
    >>> plot(x,y);
    >>> [alpha,A] = QBDQueueSTD(B, L, F, L0, true);
    >>> CheckPHRepresentation(alpha,A)
         1    
    
