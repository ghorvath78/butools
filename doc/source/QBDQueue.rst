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
  
    For Python/Numpy:
    
    >>> B=ml.matrix([[6, 1, 0], [ 0, 4, 1], [ 2, 0, 0]])
    >>> F=ml.matrix([[0, 0, 0], [ 5, 0, 0], [ 1, 3, 0]])
    >>> L=ml.matrix([[-12, 3, 2], [ 0, -14, 4], [ 3, 1, -10]])
    >>> L0 = L+B
    >>> qld,qlm = QBDQueue(B, L, F, L0, 'qlDistr', np.arange(0,11), 'qlMoms', 5)
    >>> print(qld)
    [  5.51311116e-01   2.52423571e-01   1.11087074e-01   4.81667842e-02   2.09322904e-02   9.09357429e-03   3.95071488e-03   1.71637902e-03   7.45677869e-04   3.23958385e-04   1.40743132e-04]
    >>> print(qlm)
    [0.79557486412722223, 2.0220155167474432, 7.3005182911938462, 35.016977690335295, 210.00771130457002]
    >>> alphap,Ap = QBDQueue(B, L, F, L0, 'qlDistrDPH')
    >>> print(alphap)
    [[ 0.21986817  0.12807792  0.10074279]]
    >>> print(Ap)
    [[ 0.          0.31673056  0.1241612 ]
     [ 0.          0.16525834  0.26713419]
     [ 0.          0.23369945  0.20253337]]
    >>> print(PmfFromDPH(alphap,Ap,np.arange(0,11)))
    [  5.51311116e-01   2.52423571e-01   1.11087074e-01   4.81667842e-02   2.09322904e-02   9.09357429e-03   3.95071488e-03   1.71637902e-03   7.45677869e-04   3.23958385e-04   1.40743132e-04]
    >>> print(MomentsFromDPH(alphap,Ap,5))
    [0.79557486412722234, 2.0220155167474436, 7.3005182911938489, 35.01697769033531, 210.00771130457014]
    >>> std, stm = QBDQueue(B, L, F, L0, 'stDistr', np.arange(0,1.1,0.1), 'stMoms', 5)
    >>> print(std)
    [ 0.          0.28193629  0.46807558  0.60126286  0.69980105  0.77363866  0.82922533  0.87114109  0.90276555  0.9266293   0.94463731]
    >>> print(stm)
    [0.3341414429334334, 0.23413101066642131, 0.24868918524065425, 0.35296597190069118, 0.62644746533901619]
    >>> beta, B = QBDQueue(B, L, F, L0, 'stDistrME')
    >>> print(beta)
    [[ 0.4  0.6 -0.9  0.9 -2.1  2.1 -2.1  1.2  0.9]]
    >>> print(B)
    [[-10.78786361   0.          -1.98527732   1.98527732 -15.          15.  -13.21213639  12.55037728   0.66175911]
     [  1.55762972 -13.33692411  -5.74129212   5.74129212 -10.89719372   10.89719372  -6.60606819   0.27518864   6.33087955]
     [  0.90071793  -1.94921231 -13.82752808   1.82752808  -5.62763377    3.62763377  -0.40404546   0.18345909   4.22058637]
     [  0.36644917   1.38412103  -2.82214192  -9.17785808  -5.62763377    3.62763377  -1.53644337  -3.01568095   8.55212432]
     [  1.09212509  -1.28254564  -1.3007758    5.3007758  -17.62763377   -1.37236623   2.73788071  -3.11621632   9.37833561]
     [  2.01240842  -0.98843048  -1.41106898   5.41106898  -4.4357247  -14.5642753   -1.26211929   0.88378368   9.37833561]
     [  1.58254688  -2.55348209   2.17582537   2.82417463  -2.88707654    0.88707654 -15.26211929   0.88378368   9.37833561]
     [  2.21244915  -2.4799533    2.17582537   2.82417463  -7.08909927    5.08909927  -5.26211929  -9.11621632   9.37833561]
     [  0.8238595    0.29722601   0.17582537   4.82417463  -7.65529822    5.65529822  -5.26211929   0.88378368  -0.62166439]]
    >>> print(CdfFromME(beta,B,np.arange(0,1.1,0.1)))
    [ 0.          0.28193629  0.46807558  0.60126286  0.69980105  0.77363866  0.82922533  0.87114109  0.90276555  0.9266293   0.94463731]
    >>> print(MomentsFromME(beta,B,5))
    [0.33414144293343262, 0.23413101066642034, 0.24868918524065348, 0.35296597190068874, 0.62644746533900875]

