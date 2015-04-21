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

    Returns various performane measures of a QBD queue.

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
    further parameters : 
        The rest of the function parameters specify the options
        and the performance measures to be computed.

        The supported performance measures and options in this 
        function are:
        
        +----------------+--------------------+--------------------------------------+
        | Parameter name | Input parameters   | Output                               |
        +================+====================+======================================+
        | "qlMoms"       | Number of moments  | The queue length moments             |
        +----------------+--------------------+--------------------------------------+
        | "qlDistr"      | A vector of points | The queue length distribution at     |
        |                |                    | the requested points                 |
        +----------------+--------------------+--------------------------------------+
        | "qlDistrMG"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-geometrically distributed     |
        |                |                    | queue length distribution            |
        +----------------+--------------------+--------------------------------------+
        | "qlDistrDPH"   | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-geometrically distributed     |
        |                |                    | queue length distribution, converted |
        |                |                    | to a discrete PH representation      |
        +----------------+--------------------+--------------------------------------+
        | "stMoms"       | Number of moments  | The sojourn time moments             |
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
        |                |                    | to a continuous PH representation    |
        +----------------+--------------------+--------------------------------------+
        | "prec"         | The precision      | Numerical precision to check if the  |
        |                |                    | input is valid and it is also used   |
        |                |                    | as a stopping condition when solving |
        |                |                    | the matrix-quadratic equation        |
        +----------------+--------------------+--------------------------------------+
        
        (The queue length related quantities include the customer 
        in the server, and the sojourn time related quantities 
        include the service times as well)
        
    Returns
    -------
    Ret : list of the performance measures
        Each entry of the list corresponds to a performance 
        measure requested. If there is just a single item, 
        then it is not put into a list.

    Notes
    -----
    "qlDistrMG" and "stDistrMG" behave much better numerically than 
    "qlDistrDPH" and "stDistrPH".

    Examples
    ========    
    For MATLAB:

    >>> B=[6, 1, 0; 0, 4, 1; 2, 0, 0];
    >>> F=[0, 0, 0; 5, 0, 0; 1, 3, 0];
    >>> L=[-12, 3, 2; 0, -14, 4; 3, 1, -10];
    >>> L0 = L+B;
    >>> [qld,qlm] = QBDQueue(B, L, F, L0, 'qlDistr', (0:10), 'qlMoms', 5);
    >>> qld
          0.55131      0.25242      0.11109     0.048167     0.020932    0.0090936    0.0039507    0.0017164   0.00074568   0.00032396   0.00014074
    >>> qlm
          0.79557        2.022       7.3005       35.017       210.01    
    >>> [alphap,Ap] = QBDQueue(B, L, F, L0, 'qlDistrDPH');
    >>> alphap
          0.21987      0.12808      0.10074
    >>> Ap
                0      0.31673      0.12416
                0      0.16526      0.26713
                0       0.2337      0.20253
    >>> PmfFromDPH(alphap,Ap,(0:10))'                
          0.55131      0.25242      0.11109     0.048167     0.020932    0.0090936    0.0039507    0.0017164   0.00074568   0.00032396   0.00014074
    >>> MomentsFromDPH(alphap,Ap,5)
          0.79557        2.022       7.3005       35.017       210.01
    >>> [std, stm] = QBDQueue(B, L, F, L0, 'stDistr', (0:0.1:1), 'stMoms', 5);
    >>> std
                0      0.28194      0.46808      0.60126       0.6998      0.77364      0.82923      0.87114      0.90277      0.92663      0.94464
    >>> stm
          0.33414      0.23413      0.24869      0.35297      0.62645
    >>> [beta, B] = QBDQueue(B, L, F, L0, 'stDistrME');
    >>> beta
              0.4          0.6         -0.9          0.9         -2.1          2.1         -2.1          1.2          0.9
    >>> B
          -10.788            0      -1.9853       1.9853          -15           15      -13.212        12.55      0.66176
           1.5576      -13.337      -5.7413       5.7413      -10.897       10.897      -6.6061      0.27519       6.3309
          0.90072      -1.9492      -13.828       1.8275      -5.6276       3.6276     -0.40405      0.18346       4.2206
          0.36645       1.3841      -2.8221      -9.1779      -5.6276       3.6276      -1.5364      -3.0157       8.5521
           1.0921      -1.2825      -1.3008       5.3008      -17.628      -1.3724       2.7379      -3.1162       9.3783
           2.0124     -0.98843      -1.4111       5.4111      -4.4357      -14.564      -1.2621      0.88378       9.3783
           1.5825      -2.5535       2.1758       2.8242      -2.8871      0.88708      -15.262      0.88378       9.3783
           2.2124        -2.48       2.1758       2.8242      -7.0891       5.0891      -5.2621      -9.1162       9.3783
          0.82386      0.29723      0.17583       4.8242      -7.6553       5.6553      -5.2621      0.88378     -0.62166
    >>> CdfFromME(beta,B,(0:0.1:1))'   
                0      0.28194      0.46808      0.60126       0.6998      0.77364      0.82923      0.87114      0.90277      0.92663      0.94464
    >>> MomentsFromME(beta,B,5)
          0.33414      0.23413      0.24869      0.35297      0.62645
  

