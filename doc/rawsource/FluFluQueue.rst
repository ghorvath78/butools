butools.queues.FluFluQueue
==========================

.. currentmodule:: butools.queues

.. np:function:: FluFluQueue

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Ret = FluFluQueue(Qin, Rin, Qout, Rout, srv0stop, ...)`
        * - Mathematica:
          - :code:`Ret = FluFluQueue[Qin, Rin, Qout, Rout, srv0stop, ...]`
        * - Python/Numpy:
          - :code:`Ret = FluFluQueue(Qin, Rin, Qout, Rout, srv0stop, ...)`

    Returns various performane measures of a fluid queue
    with independent fluid arrival and service processes.

    Two types of boundary behavior is available. If 
    srv0stop=false, the output process evolves continuously
    even if the queue is empty. If srv0stop=true, the 
    output process slows down if there is fewer fluid in
    the queue than it can serve. If the queue is empty
    and the fluid input rate is zero, the output process
    freezes till fluid arrives.

    Parameters
    ----------
    Qin : matrix, shape (N,N)
        The generator of the background Markov chain 
        corresponding to the input process
    Rin : matrix, shape (N,N)
        Diagonal matrix containing the fluid input rates
        associated to the states of the input background 
        process
    Qout : matrix, shape (N,N)
        The generator of the background Markov chain 
        corresponding to the output process
    Rout : matrix, shape (N,N)
        Diagonal matrix containing the fluid output rates
        associated to the states of the input background 
        process
    srv0stop : bool
        If true, the service output process slows down if
        there is fewer fluid in the queue than it can 
        serve. If false, the output process evolves 
        continuously.
    further parameters : 
        The rest of the function parameters specify the options
        and the performance measures to be computed.

        The supported performance measures and options in this 
        function are:
        
        +----------------+--------------------+--------------------------------------+
        | Parameter name | Input parameters   | Output                               |
        +================+====================+======================================+
        | "flMoms"       | Number of moments  | The moments of the fluid level       |
        +----------------+--------------------+--------------------------------------+
        | "flDistr"      | A vector of points | The fluid level distribution at      |
        |                |                    | the requested points (cdf)           |
        +----------------+--------------------+--------------------------------------+
        | "flDistrME"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-exponentially distributed     |
        |                |                    | fluid level distribution             |
        +----------------+--------------------+--------------------------------------+
        | "flDistrPH"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-exponentially distributed     |
        |                |                    | fluid level distribution, converted  |
        |                |                    | to a PH representation               |
        +----------------+--------------------+--------------------------------------+
        | "stMoms"       | Number of moments  | The sojourn time moments of fluid    |
        |                |                    | drops                                |
        +----------------+--------------------+--------------------------------------+
        | "stDistr"      | A vector of points | The sojourn time distribution at the |
        |                |                    | requested points (cummulative, cdf)  |
        +----------------+--------------------+--------------------------------------+
        | "stDistrME"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-exponentially distributed     |
        |                |                    | sojourn time distribution            |
        +----------------+--------------------+--------------------------------------+
        | "stDistrPH"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-exponentially distributed     |
        |                |                    | sojourn time distribution, converted |
        |                |                    | to a PH representation               |
        +----------------+--------------------+--------------------------------------+
        | "prec"         | The precision      | Numerical precision to check if the  |
        |                |                    | input is valid and it is also used   |
        |                |                    | as a stopping condition when solving |
        |                |                    | the Riccati equation                 |
        +----------------+--------------------+--------------------------------------+

    Returns
    -------
    Ret : list of the performance measures
        Each entry of the list corresponds to a performance 
        measure requested. If there is just a single item, 
        then it is not put into a list.

    Notes
    -----
    "flDistrME" and "stDistrME" behave much better numerically than 
    "flDistrPH" and "stDistrPH".

    References
    ----------
    .. [1] Horvath G, Telek M, "Sojourn times in fluid queues 
           with independent and dependent input and output 
           processes PERFORMANCE EVALUATION 79: pp. 160-181, 2014.

    Examples
    ========    

