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
    "qlDistrME" and "stDistrME" behave much better numerically than 
    "qlDistrPH" and "stDistrPH".
    
    Examples
    ========    
    For MATLAB:

    >>> Q = [-9 2 4 0 1 2; 6 -25 5 3 7 4; 1 3 -4 0 0 0; 0 0 0 -8 3 5; 7 3 0 2 -13 1; 7 8 0 3 8 -26];
    >>> Rin = diag([4 2 1 0 0 3]);
    >>> Rout = diag([6 2 0 0 3 2]);
    >>> [qld,qlm] = FluidQueue(Q, Rin, Rout, 'qlDistr', (0:0.5:2), 'qlMoms', 5);
    >>> qld
          0.23662      0.70447      0.87861      0.95014      0.97952
    >>> qlm
          0.40636       0.4546      0.76609        1.722       4.8384
    >>> [alphap,Ap] = FluidQueue(Q, Rin, Rout, 'qlDistrPH');
    >>> alphap
          0.63124      0.13213
    >>> Ap
          -2.0387      0.41483
             12.1      -21.143
    >>> CdfFromPH(alphap,Ap,(0:0.5:2))'                
          0.23662      0.70447      0.87861      0.95014      0.97952
    >>> MomentsFromPH(alphap,Ap,5)
          0.40636       0.4546      0.76609        1.722       4.8384
    >>> [std, stm] = FluidQueue(Q, Rin, Rout, 'stDistr', (0:0.2:1), 'stMoms', 5);
    >>> std
          0.31678      0.67546       0.7949       0.8688      0.91618      0.94652
    >>> stm
          0.23252      0.20069      0.26684      0.47402       1.0523
    >>> [beta, B] = FluidQueue(Q, Rin, Rout, 'stDistrME');
    >>> beta
       6.1487e-17            0     0.078831    -0.087519     0.095339   1.1102e-16            0     -0.10977      0.10977     -0.10977      0.10977      0.59657
    >>> B
          -22.232      -9.0243       12.676      -16.419      -10.249       30.265            0            0            0      -20.594       20.594            0
          -7.8856          -14       12.876       -9.097      -5.6783       16.769            0            0            0       -11.41        11.41            0
          -7.2098       15.443      -18.077      -16.653       4.8244       14.325      -4.9559       5.1109            0      -13.925       13.925            0
          -6.4942        13.91       45.225      -88.757       8.1259       15.247      -9.3289       9.6206      -6.2713      -6.2713       11.288       1.2543
          -5.1787       9.9439       41.788      -50.497      -22.618       13.996      -8.5637       3.0746            0      -6.9083       11.514       1.1514
           1.1565       9.2439       33.466      -40.969       105.27      -123.85       -13.15       9.2801            0       -6.422        3.211       8.5627
           1.1143       8.9065       32.244      -39.473       97.844      -111.62       -21.67       14.098            0      -6.1876       3.0938       8.2501
           1.0805       5.2804        29.45      -31.898        93.14      -109.17       47.142      -53.616            0           -6            0           11
           1.0805       3.4399       31.604       -34.29       95.746      -109.17       47.142      -45.616           -8           -8            2           11
           1.0805       3.4399        29.45      -31.898       95.746      -109.17       44.233      -42.616            0          -16            0           13
           1.4204     -0.85454       34.477      -38.276       102.69      -109.17       44.233      -45.616            3       3.1162      -24.622       18.505
           1.4204     -0.85454       28.732      -31.898       101.83      -108.24       37.445      -38.616            0      -2.0893       60.013      -57.924
    >>> CdfFromME(beta,B,(0:0.2:1))'   
          0.31678      0.67546       0.7949       0.8688      0.91618      0.94652
    >>> MomentsFromME(beta,B,5)
          0.23252      0.20069      0.26684      0.47402       1.0523
 
