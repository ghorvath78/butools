butools.queues.QBDQueue
=======================

.. currentmodule:: butools.queues

.. np:function:: QBDQueue

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Ret = QBDQueue(B, L, F, L0, ...)`
        * - Mathematica:
          - :code:`Ret = QBDQueue[B, L, F, L0, ...]`
        * - Python/Numpy:
          - :code:`Ret = QBDQueue(B, L, F, L0, ...)`

    Returns various performane measures of a continuous time
    QBD queue.

    QBD queues have a background continuous time Markov chain
    with generator Q whose the transitions can be partitioned
    into three sets: transitions accompanied by an arrival
    of a new job (F, forward), transitions accompanied by 
    the service of the current job in the server (B, 
    backward) and internal transitions (L, local). 
    Thus we have Q=B+L+F. L0 is the matrix of local 
    transition rates if the queue is empty.

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
    further parameters : 
        The rest of the function parameters specify the options
        and the performance measures to be computed.

        The supported performance measures and options in this 
        function are:
        
        +----------------+--------------------+----------------------------------------+
        | Parameter name | Input parameters   | Output                                 |
        +================+====================+========================================+
        | "ncMoms"       | Number of moments  | The moments of the number of customers |
        +----------------+--------------------+----------------------------------------+
        | "ncDistr"      | Upper limit K      | The distribution of the number of      |
        |                |                    | customers from level 0 to level K-1    |
        +----------------+--------------------+----------------------------------------+
        | "ncDistrMG"    | None               | The vector-matrix parameters of the    |
        |                |                    | matrix-geometric distribution of the   |
        |                |                    | number of customers in the system      |
        +----------------+--------------------+----------------------------------------+
        | "ncDistrDPH"   | None               | The vector-matrix parameters of the    |
        |                |                    | matrix-geometric distribution of the   |
        |                |                    | number of customers in the system,     |
        |                |                    | converted to a discrete PH             |
        |                |                    | representation                         |
        +----------------+--------------------+----------------------------------------+
        | "stMoms"       | Number of moments  | The sojourn time moments               |
        +----------------+--------------------+----------------------------------------+
        | "stDistr"      | A vector of points | The sojourn time distribution at the   |
        |                |                    | requested points (cummulative, cdf)    |
        +----------------+--------------------+----------------------------------------+
        | "stDistrME"    | None               | The vector-matrix parameters of the    |
        |                |                    | matrix-exponentially distributed       |
        |                |                    | sojourn time distribution              |
        +----------------+--------------------+----------------------------------------+
        | "stDistrPH"    | None               | The vector-matrix parameters of the    |
        |                |                    | matrix-exponentially distributed       |
        |                |                    | sojourn time distribution, converted   |
        |                |                    | to a continuous PH representation      |
        +----------------+--------------------+----------------------------------------+
        | "prec"         | The precision      | Numerical precision used as a stopping |
        |                |                    | condition when solving the             |
        |                |                    | matrix-quadratic equation              |
        +----------------+--------------------+----------------------------------------+
        
        (The quantities related to the number of customers in 
        the system include the customer in the server, and the 
        sojourn time related quantities include the service 
        times as well)
        
    Returns
    -------
    Ret : list of the performance measures
        Each entry of the list corresponds to a performance 
        measure requested. If there is just a single item, 
        then it is not put into a list.

    Notes
    -----
    "ncDistrMG" and "stDistrMG" behave much better numerically than 
    "ncDistrDPH" and "stDistrPH".

    Examples
    ========    
    
