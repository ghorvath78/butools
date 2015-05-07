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
        | "qlMoms"       | Number of moments  | The moments of the fluid level       |
        +----------------+--------------------+--------------------------------------+
        | "qlDistr"      | A vector of points | The fluid level distribution at      |
        |                |                    | the requested points (cdf)           |
        +----------------+--------------------+--------------------------------------+
        | "qlDistrME"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-exponentially distributed     |
        |                |                    | fluid level distribution             |
        +----------------+--------------------+--------------------------------------+
        | "qlDistrPH"    | None               | The vector-matrix parameters of the  |
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
    "qlDistrME" and "stDistrME" behave much better numerically than 
    "qlDistrPH" and "stDistrPH".

    Examples
    ========    
    For MATLAB:
    
    >>> Qin = [-2 1 1; 2 -5 3; 4 0 -4];
    >>> Rin = diag([3 7 0]);
    >>> Qout = [-4 1 3; 6 -8 2; 3 7 -10];
    >>> Rout = diag([1 7 15]);
    >>> [qld, qlm, std] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'qlDistr', (0:0.5:2), 'qlMoms', 5, 'stDistrPH');
    >>> qld
           0.3918      0.68193      0.81605      0.88805      0.93027
    >>> qlm
           0.5357       1.0765       3.4298        14.87       81.162
    >>> [beta, B] = std{:};
    >>> beta
          0.45285      0.24877
    >>> B
          -5.4973      0.83675
           1.4492      -3.7914
    >>> [std, qld, stm] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'stDistr', (0:0.5:2), 'qlDistrME', 'stMoms', 5);
    >>> std
          0.57864      0.92608      0.98589      0.99725      0.99946
    >>> stm
          0.12096     0.071546     0.064592      0.07852      0.11997
    >>> [alpha, A] = qld{:};
    >>> alpha
         -0.24261      0.59524
    >>> A
          -2.1425       1.5194
          0.43807      -1.4538

    For Python/Numpy:
    
    >>> Qin = ml.matrix([[-2, 1, 1], [ 2, -5, 3], [ 4, 0, -4]])
    >>> Rin = Diag(ml.matrix([[3, 7, 0]]))
    >>> Qout = ml.matrix([[-4, 1, 3], [ 6, -8, 2], [ 3, 7, -10]])
    >>> Rout = Diag(ml.matrix([[1, 7, 15]]))
    >>> qld, qlm, std = FluFluQueue(Qin, Rin, Qout, Rout, False, 'qlDistr', np.arange(0,2.5,0.5), 'qlMoms', 5, 'stDistrPH')
    >>> print(qld)
    [ 0.39180087  0.68192691  0.81605281  0.88804626  0.93026677]
    >>> print(qlm)
    [0.53570356276760078, 1.0765385576008903, 3.4298090557057193, 14.869885651621992, 81.162335691323889]
    >>> beta, B = std
    >>> print(beta)
    [[ 0.45285051  0.24877407]]
    >>> print(B)
    [[-5.49734397  0.8367526 ]
     [ 1.44916617 -3.79142652]]
    >>> std, qld, stm = FluFluQueue(Qin, Rin, Qout, Rout, True, 'stDistr', np.arange(0,2.5,0.5), 'qlDistrME', 'stMoms', 5)
    >>> print(std)
    [[ 0.57863816  0.92608217  0.985895    0.99724611  0.99945918]]
    >>> print(stm)
    [0.12096271585043798, 0.071546261021158727, 0.064592069915750366, 0.078520337012782174, 0.11996554734440762]
    >>> alpha, A = qld
    >>> print(alpha)
    [[-0.24260694  0.59524366]]
    >>> print(A)
    [[-2.14252917  1.51937881]
     [ 0.4380681  -1.45384643]]

