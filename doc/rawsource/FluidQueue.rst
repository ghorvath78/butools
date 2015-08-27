butools.queues.FluidQueue
=========================

.. currentmodule:: butools.queues

.. np:function:: FluidQueue

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Ret = FluidQueue(Q, Rin, Rout, ...)`
        * - Mathematica:
          - :code:`Ret = FluidQueue[Q, Rin, Rout, ...]`
        * - Python/Numpy:
          - :code:`Ret = FluidQueue(Q, Rin, Rout, ...)`

    Returns various performane measures of a fluid queue.

    In a fluid queue there is a background continuous time
    Markov chain (given by generator Q), and diagonal
    matrix Rin (Rout) whose ith entry provides the 
    fluid rate at which fluid enters the queue (can be 
    served) while the background process is in state i.

    Parameters
    ----------
    Q : matrix, shape (N,N)
        The generator of the background Markov chain
    Rin : matrix, shape (N,N)
        Diagonal matrix containing the fluid input rates
        associated to the states of the background process
    Rout : matrix, shape (N,N)
        Diagonal matrix containing the fluid output rates
        associated to the states of the background process
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
        | "Q0"           | Matrix, shape(N,N) | The generator of the background      |
        |                |                    | Markov chain when the fluid level is |
        |                |                    | zero. If not given, Q0=Q is assumed  |
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
    
    Examples
    ========    

