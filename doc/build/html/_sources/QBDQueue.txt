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
    For Matlab:

    >>> B = [6., 1., 0.; 0., 4., 1.; 2., 0., 0.];
    >>> F = [0., 1., 1.; 5., 0., 0.; 1., 3., 0.];
    >>> L = [-14., 3., 2.; 0., -14., 4.; 3., 1., -10.];
    >>> L0 = L+B;
    >>> disp(L0);
        -8     4     2
         0   -10     5
         5     1   -10
    >>> [qld, qlm] = QBDQueue(B, L, F, L0, 'qlDistr', (0:1:10), 'qlMoms', 5);
    >>> disp(qld);
      Columns 1 through 6
          0.29094      0.20433      0.14547      0.10358     0.073709     0.052461
      Columns 7 through 11
         0.037338     0.026574     0.018913     0.013461    0.0095805
    >>> disp(qlm);
             2.46       14.608       128.88       1515.9        22288
    >>> [alphap, Ap] = QBDQueue(B, L, F, L0, 'qlDistrDPH');
    >>> disp(alphap);
          0.28256      0.22386      0.20264
    >>> disp(Ap);
          0.11832       0.4103      0.18373
            0.196      0.19477      0.32261
          0.27566      0.22146      0.21229
    >>> [alpha, A] = QBDQueue(B, L, F, L0, 'qlDistrMG');
    >>> disp(alpha);
         0.042581     -0.27854      0.94503
    >>> disp(A);
       -0.0056928     -0.63242       1.3313
         0.069475     -0.36133      0.99875
         0.037353     -0.21854       0.8924
    >>> qldFromDPH = PmfFromDPH(alphap, Ap, (0:1:10));
    >>> disp(qldFromDPH);
      Columns 1 through 6
          0.29094      0.20433      0.14547      0.10358     0.073709     0.052461
      Columns 7 through 11
         0.037338     0.026574     0.018913     0.013461    0.0095805
    >>> qlmFromMG = MomentsFromMG(alpha, A, 5);
    >>> disp(qlmFromMG);
             2.46       14.608       128.88       1515.9        22288
    >>> [std, stm] = QBDQueue(B, L, F, L0, 'stDistr', (0.:0.1:1.), 'stMoms', 5);
    >>> disp(std);
      Columns 1 through 6
       1.1102e-16      0.14225      0.25715      0.35489      0.43933      0.51262
      Columns 7 through 11
           0.5763      0.63165      0.67977      0.72161      0.75797
    >>> disp(stm);
          0.70236       1.0017       2.1463       6.1325       21.903
    >>> [betap, Bp] = QBDQueue(B, L, F, L0, 'stDistrPH');
    >>> disp(betap);
      Columns 1 through 6
          0.52899            0            0            0      0.35507            0
      Columns 7 through 9
                0            0      0.11594
    >>> disp(Bp);
      Columns 1 through 6
          -12.591        2.362      0.50976       4.2349      0.39366     0.084961
            2.001      -12.662      0.95667       0.3335       4.2231      0.15944
           2.8354      0.93483      -13.595      0.47257      0.15581       4.0674
                5            0            0      -13.061       1.5746      0.33984
                0            5            0        1.334      -13.108      0.63778
                0            0            5       1.8903      0.62322       -13.73
           4.4697      0.78732      0.16992            4            0            0
            0.667       4.4461      0.31889            0            4            0
          0.94514      0.31161       4.1349            0            0            4
      Columns 7 through 9
                3            0            0
                0            3            0
                0            0            3
           4.2349      0.39366     0.084961
           0.3335       4.2231      0.15944
          0.47257      0.15581       4.0674
              -10            0            0
                0          -10            0
                0            0          -10
    >>> [beta, B] = QBDQueue(B, L, F, L0, 'stDistrME');
    >>> disp(beta);
      Columns 1 through 6
          0.17391      0.47826     -0.71739      0.71739      -1.2391       1.2391
      Columns 7 through 9
          -1.2391      0.52174       1.0652
    >>> disp(B);
      Columns 1 through 6
          -12.591     -0.41431       -3.408        2.165      -13.136           15
           0.5517      -13.315      -8.0386       7.3408       -10.01       10.714
          0.20504      -1.8844      -15.359       1.8939      -4.0612       3.4626
         -0.28868       1.3811      -2.5957      -11.073      -3.7559       3.4626
          0.48659      -1.3755      -1.3906       5.4522      -17.351      -1.5374
           1.5997     -0.77852      -1.6145       5.5725      -4.1872      -15.012
           1.0779      -2.1431       1.8972       1.9933      -1.6703      0.26853
           2.0033      -1.9938       1.8972       1.9933      -5.9052       4.3998
          0.59557      0.82156      -1.1028       4.9933      -6.4916       4.9354
      Columns 7 through 9
          -13.409       12.066       1.3431
          -6.3994     0.033031       6.6716
         -0.19879     0.022021       4.4477
          -1.3717      -3.2512       8.8939
            2.853      -3.4525       9.8705
          -0.7327       0.5475       9.8705
          -14.463       0.5475       9.8705
          -4.3594      -9.4525       9.8705
          -4.3085       0.5475     -0.12955
    >>> stdFromPH = CdfFromPH(betap, Bp, (0.:0.1:1.));
    >>> disp(stdFromPH);
      Columns 1 through 6
       1.1102e-16      0.14225      0.25715      0.35489      0.43933      0.51262
      Columns 7 through 11
           0.5763      0.63165      0.67977      0.72161      0.75797
    >>> stmFromME = MomentsFromME(beta, B, 5);
    >>> disp(stmFromME);
          0.70236       1.0017       2.1463       6.1325       21.903

    For Mathematica:

    
    For Python/Numpy:

    >>> B = ml.matrix([[6., 1., 0.],[0., 4., 1.],[2., 0., 0.]])
    >>> F = ml.matrix([[0., 1., 1.],[5., 0., 0.],[1., 3., 0.]])
    >>> L = ml.matrix([[-14., 3., 2.],[0., -14., 4.],[3., 1., -10.]])
    >>> L0 = L+B
    >>> print(L0)
    [[ -8.   4.   2.]
     [  0. -10.   5.]
     [  5.   1. -10.]]
    >>> qld, qlm = QBDQueue(B, L, F, L0, "qlDistr", np.arange(0,11.0,1), "qlMoms", 5)
    Final Residual Error for G:  9.71445146547e-17
    Final Residual Error for R:  1.11022302463e-16
    >>> print(qld)
    [ 0.29094  0.20433  0.14547  0.10358  0.07371  0.05246  0.03734  0.02657  0.01891  0.01346  0.00958]
    >>> print(qlm)
    [2.4600295202120317, 14.607919051794482, 128.87860443055746, 1515.8884540468425, 22287.969412532377]
    >>> alphap, Ap = QBDQueue(B, L, F, L0, "qlDistrDPH")
    Final Residual Error for G:  9.71445146547e-17
    Final Residual Error for R:  1.11022302463e-16
    >>> print(alphap)
    [[ 0.28256  0.22386  0.20264]]
    >>> print(Ap)
    [[ 0.11832  0.4103   0.18373]
     [ 0.196    0.19477  0.32261]
     [ 0.27566  0.22146  0.21229]]
    >>> alpha, A = QBDQueue(B, L, F, L0, "qlDistrMG")
    Final Residual Error for G:  9.71445146547e-17
    Final Residual Error for R:  1.11022302463e-16
    >>> print(alpha)
    [[ 0.04258 -0.27854  0.94503]]
    >>> print(A)
    [[-0.00569 -0.63242  1.33129]
     [ 0.06947 -0.36133  0.99875]
     [ 0.03735 -0.21854  0.8924 ]]
    >>> qldFromDPH = PmfFromDPH(alphap, Ap, np.arange(0,11.0,1))
    >>> print(qldFromDPH)
    [ 0.29094  0.20433  0.14547  0.10358  0.07371  0.05246  0.03734  0.02657  0.01891  0.01346  0.00958]
    >>> qlmFromMG = MomentsFromMG(alpha, A, 5)
    >>> print(qlmFromMG)
    [2.4600295202120313, 14.607919051794475, 128.87860443055737, 1515.8884540468407, 22287.969412532337]
    >>> std, stm = QBDQueue(B, L, F, L0, "stDistr", np.arange(0.,1.1,0.1), "stMoms", 5)
    Final Residual Error for G:  9.71445146547e-17
    Final Residual Error for R:  1.11022302463e-16
    >>> print(std)
    [  1.11022e-16   1.42250e-01   2.57146e-01   3.54886e-01   4.39333e-01   5.12619e-01   5.76300e-01   6.31653e-01   6.79773e-01   7.21607e-01   7.57975e-01]
    >>> print(stm)
    [0.70235625432140525, 1.0017175985977307, 2.1462963079180937, 6.1325315761778585, 21.903161812239034]
    >>> betap, Bp = QBDQueue(B, L, F, L0, "stDistrPH")
    Final Residual Error for G:  9.71445146547e-17
    Final Residual Error for R:  1.11022302463e-16
    >>> print(betap)
    [[ 0.52899  0.       0.       0.       0.35507  0.       0.       0.       0.11594]]
    >>> print(Bp)
    [[-12.59079   2.36197   0.50976   4.23487   0.39366   0.08496   3.        0.        0.     ]
     [  2.00101 -12.66159   0.95667   0.3335    4.22307   0.15944   0.        3.        0.     ]
     [  2.83542   0.93483 -13.59536   0.47257   0.15581   4.06744   0.        0.        3.     ]
     [  5.        0.        0.      -13.06053   1.57465   0.33984   4.23487   0.39366   0.08496]
     [  0.        5.        0.        1.33401 -13.10772   0.63778   0.3335    4.22307   0.15944]
     [  0.        0.        5.        1.89028   0.62322 -13.73024   0.47257   0.15581   4.06744]
     [  4.46974   0.78732   0.16992   4.        0.        0.      -10.        0.        0.     ]
     [  0.667     4.44614   0.31889   0.        4.        0.        0.      -10.        0.     ]
     [  0.94514   0.31161   4.13488   0.        0.        4.        0.        0.      -10.     ]]
    >>> beta, B = QBDQueue(B, L, F, L0, "stDistrME")
    Final Residual Error for G:  9.71445146547e-17
    Final Residual Error for R:  1.11022302463e-16
    >>> print(beta)
    [[ 0.17391  0.47826 -0.71739  0.71739 -1.23913  1.23913 -1.23913  0.52174  1.06522]]
    >>> print(B)
    [[-12.59079  -0.41431  -3.40797   2.16504 -13.13561  15.      -13.40921  12.06606   1.34314]
     [  0.5517  -13.31488  -8.03859   7.34081 -10.01035  10.7136   -6.39935   0.03303   6.67157]
     [  0.20504  -1.88439 -15.35906   1.89388  -4.0612    3.46259  -0.19879   0.02202   4.44771]
     [ -0.28868   1.38111  -2.59572 -11.07296  -3.75595   3.46259  -1.37175  -3.25117   8.89385]
     [  0.48659  -1.37547  -1.39064   5.45219 -17.35131  -1.53741   2.85299  -3.4525    9.87045]
     [  1.59973  -0.77852  -1.61449   5.57247  -4.1872  -15.01224  -0.7327    0.5475    9.87045]
     [  1.07788  -2.14306   1.89721   1.99333  -1.67029   0.26853 -14.46294   0.5475    9.87045]
     [  2.00326  -1.99382   1.89721   1.99333  -5.90516   4.39982  -4.35937  -9.4525    9.87045]
     [  0.59557   0.82156  -1.10279   4.99333  -6.49163   4.93542  -4.30849   0.5475   -0.12955]]
    >>> stdFromPH = CdfFromPH(betap, Bp, np.arange(0.,1.1,0.1))
    >>> print(stdFromPH)
    [  1.11022e-16   1.42250e-01   2.57146e-01   3.54886e-01   4.39333e-01   5.12619e-01   5.76300e-01   6.31653e-01   6.79773e-01   7.21607e-01   7.57975e-01]
    >>> stmFromME = MomentsFromME(beta, B, 5)
    >>> print(stmFromME)
    [0.70235625432140614, 1.001717598597732, 2.1462963079180968, 6.1325315761778718, 21.903161812239095]

