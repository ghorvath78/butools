butools.queues.MMAPPH1FCFS
==========================

.. currentmodule:: butools.queues

.. np:function:: MMAPPH1FCFS

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Ret = MMAPPH1FCFS(D, sigma, S, ...)`
        * - Mathematica:
          - :code:`Ret = MMAPPH1FCFS[D, sigma, S, ...]`
        * - Python/Numpy:
          - :code:`Ret = MMAPPH1FCFS(D, sigma, S, ...)`

    Returns various performane measures of a MMAP[K]/PH[K]/1 
    first-come-first-serve queue, see [1]_.

    Parameters
    ----------
    D : list of matrices of shape (N,N), length (K+1)
        The D0...DK matrices of the arrival process.
    sigma : list of row vectors, length (K)
        The list containing the initial probability vectors of the service
        time distributions of the various customer types. The length of the
       vectors does not have to be the same.
    S : list of square matrices, length (K)
        The transient generators of the phase type distributions representing
        the service time of the jobs belonging to various types.
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
        |                |                    | condition when solving the Riccati     |
        |                |                    | equation                               |
        +----------------+--------------------+----------------------------------------+
        | "classes"      | Vector of integers | Only the performance measures          |
        |                |                    | belonging to these classes are         |
        |                |                    | returned. If not given, all classes    |
        |                |                    | are analyzed.                          |
        +----------------+--------------------+----------------------------------------+
        
        (The quantities related to the number of customers in 
        the system include the customer in the server, and the 
        sojourn time related quantities include the service 
        times as well)

    Returns
    -------
    Ret : list of the performance measures
        Each entry of the list corresponds to a performance 
        measure requested. Each entry is a matrix, where the
        columns belong to the various job types.
        If there is just a single item, 
        then it is not put into a list.

    References
    ----------
    .. [1] Qiming He, "Analysis of a continuous time 
           SM[K]/PH[K]/1/FCFS queue: Age process, sojourn times,
           and queue lengths", Journal of Systems Science and 
           Complexity, 25(1), pp 133-155, 2012.

    Examples
    ========    

