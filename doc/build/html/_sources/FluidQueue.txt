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
    For Matlab:

    >>> Q = [-9., 2., 4., 0., 1., 2.; 6., -25., 5., 3., 7., 4.; 1., 3., -4., 0., 0., 0.; 0., 0., 0., -8., 3., 5.; 7., 3., 0., 2., -13., 1.; 7., 8., 0., 3., 8., -26.];
    >>> vRin = [4.,2.,1.,0.,0.,3.];
    >>> vRout = [6.,2.,0.,0.,3.,2.];
    >>> Rin = diag(vRin);
    >>> Rout = diag(vRout);
    >>> [qld, qlm] = FluidQueue(Q, Rin, Rout, 'qlDistr', (0.:0.1:1.), 'qlMoms', 5);
    Final Residual Error for Psi:    5.3291e-15
    >>> disp(qld);
      Columns 1 through 6
          0.23662      0.39265      0.49537      0.57807       0.6469      0.70447
      Columns 7 through 11
          0.75265      0.79297      0.82672      0.85497      0.87861
    >>> disp(qlm);
          0.40636       0.4546      0.76609        1.722       4.8384
    >>> [alphap, Ap] = FluidQueue(Q, Rin, Rout, 'qlDistrPH');
    Final Residual Error for Psi:    5.3291e-15
    >>> disp(alphap);
          0.63124      0.13213
    >>> disp(Ap);
          -2.0387      0.41483
             12.1      -21.143
    >>> [alpha, A] = FluidQueue(Q, Rin, Rout, 'qlDistrME');
    Final Residual Error for Psi:    5.3291e-15
    >>> disp(alpha);
         0.099033      0.66434
    >>> disp(A);
          -3.8739       2.2653
           16.206      -19.308
    >>> qldFromPH = CdfFromPH(alphap, Ap, (0.:0.1:1.));
    >>> disp(qldFromPH);
      Columns 1 through 6
          0.23662      0.39265      0.49537      0.57807       0.6469      0.70447
      Columns 7 through 11
          0.75265      0.79297      0.82672      0.85497      0.87861
    >>> qlmFromME = MomentsFromME(alpha, A, 5);
    >>> disp(qlmFromME);
          0.40636       0.4546      0.76609        1.722       4.8384
    >>> [std, stm] = FluidQueue(Q, Rin, Rout, 'stDistr', (0.:0.1:1.), 'stMoms', 5);
    Final Residual Error for Psi:    5.3291e-15
    >>> disp(std);
      Columns 1 through 6
          0.31678      0.57513      0.67546      0.74313       0.7949      0.83595
      Columns 7 through 11
           0.8688      0.89511      0.91618      0.93304      0.94652
    >>> disp(stm);
          0.23252      0.20069      0.26684      0.47402       1.0523
    >>> [betap, Bp] = FluidQueue(Q, Rin, Rout, 'stDistrPH');
    Final Residual Error for Psi:    5.3291e-15
    >>> disp(betap);
      Columns 1 through 6
          0.24893     0.030403     0.068594     0.023408      0.21436            0
      Columns 7 through 12
                0            0            0            0            0     0.097532
    >>> disp(Bp);
      Columns 1 through 6
          -21.232        2.489            2            0            4            0
             72.6      -135.86            0            2            0            4
                6            0      -29.077      0.82967            5            0
                0            6         24.2      -67.286            0            5
                1            0            3            0           -4            0
                0            1            0            3            0           -4
                0            0            0            0            0            0
                0            0            0            0            0            0
                7            0            3            0            0            0
                0            7            0            3            0            0
                7            0            8            0            0            0
                0            7            0            8            0            0
      Columns 7 through 12
                0            0            1            0            2            0
                0            0            0            1            0            2
                3            0            7            0            4            0
                0            3            0            7            0            4
                0            0            0            0            0            0
                0            0            0            0            0            0
               -8            0            3            0            5            0
                0           -8            0            3            0            5
                2            0      -19.116       1.2445            1            0
                0            2         36.3      -76.429            0            1
                3            0            8            0      -30.077      0.82967
                0            3            0            8         24.2      -68.286
    >>> [beta, B] = FluidQueue(Q, Rin, Rout, 'stDistrME');
    Final Residual Error for Psi:    5.3291e-15
    >>> disp(beta);
      Columns 1 through 6
       6.1487e-17            0     0.078831    -0.087519     0.095339   1.1102e-16
      Columns 7 through 12
                0     -0.10977      0.10977     -0.10977      0.10977      0.59657
    >>> disp(B);
      Columns 1 through 6
          -22.232      -9.0243       12.676      -16.419      -10.249       30.265
          -7.8856          -14       12.876       -9.097      -5.6783       16.769
          -7.2098       15.443      -18.077      -16.653       4.8244       14.325
          -6.4942        13.91       45.225      -88.757       8.1259       15.247
          -5.1787       9.9439       41.788      -50.497      -22.618       13.996
           1.1565       9.2439       33.466      -40.969       105.27      -123.85
           1.1143       8.9065       32.244      -39.473       97.844      -111.62
           1.0805       5.2804        29.45      -31.898        93.14      -109.17
           1.0805       3.4399       31.604       -34.29       95.746      -109.17
           1.0805       3.4399        29.45      -31.898       95.746      -109.17
           1.4204     -0.85454       34.477      -38.276       102.69      -109.17
           1.4204     -0.85454       28.732      -31.898       101.83      -108.24
      Columns 7 through 12
                0            0            0      -20.594       20.594            0
                0            0            0       -11.41        11.41            0
          -4.9559       5.1109            0      -13.925       13.925            0
          -9.3289       9.6206      -6.2713      -6.2713       11.288       1.2543
          -8.5637       3.0746            0      -6.9083       11.514       1.1514
           -13.15       9.2801            0       -6.422        3.211       8.5627
           -21.67       14.098            0      -6.1876       3.0938       8.2501
           47.142      -53.616            0           -6            0           11
           47.142      -45.616           -8           -8            2           11
           44.233      -42.616            0          -16            0           13
           44.233      -45.616            3       3.1162      -24.622       18.505
           37.445      -38.616            0      -2.0893       60.013      -57.924
    >>> stdFromPH = CdfFromPH(betap, Bp, (0.:0.1:1.));
    >>> disp(stdFromPH);
      Columns 1 through 6
          0.31678      0.57513      0.67546      0.74313       0.7949      0.83595
      Columns 7 through 11
           0.8688      0.89511      0.91618      0.93304      0.94652
    >>> stmFromME = MomentsFromME(beta, B, 5);
    >>> disp(stmFromME);
          0.23252      0.20069      0.26684      0.47402       1.0523

    For Mathematica:

    
    For Python/Numpy:

    >>> Q = ml.matrix([[-9., 2., 4., 0., 1., 2.],[6., -25., 5., 3., 7., 4.],[1., 3., -4., 0., 0., 0.],[0., 0., 0., -8., 3., 5.],[7., 3., 0., 2., -13., 1.],[7., 8., 0., 3., 8., -26.]])
    >>> vRin = ml.matrix([[4.,2.,1.,0.,0.,3.]])
    >>> vRout = ml.matrix([[6.,2.,0.,0.,3.,2.]])
    >>> Rin = Diag(vRin)
    >>> Rout = Diag(vRout)
    >>> qld, qlm = FluidQueue(Q, Rin, Rout, "qlDistr", np.arange(0.,1.1,0.1), "qlMoms", 5)
    Final Residual Error for G:  6.661338147750939e-15
    >>> print(qld)
    [ 0.23662  0.39265  0.49537  0.57807  0.6469   0.70447  0.75265  0.79297  0.82672  0.85497  0.87861]
    >>> print(qlm)
    [0.406361527249128, 0.45459875335830258, 0.76609389063916067, 1.7219814531396287, 4.8383553122056897]
    >>> alphap, Ap = FluidQueue(Q, Rin, Rout, "qlDistrPH")
    Final Residual Error for G:  6.661338147750939e-15
    >>> print(alphap)
    [[ 0.63124  0.13213]]
    >>> print(Ap)
    [[ -2.03873   0.41483]
     [ 12.10002 -21.1431 ]]
    >>> alpha, A = FluidQueue(Q, Rin, Rout, "qlDistrME")
    Final Residual Error for G:  6.661338147750939e-15
    >>> print(alpha)
    [[ 0.09903  0.66434]]
    >>> print(A)
    [[ -3.87389   2.26527]
     [ 16.20615 -19.30794]]
    >>> qldFromPH = CdfFromPH(alphap, Ap, np.arange(0.,1.1,0.1))
    >>> print(qldFromPH)
    [ 0.23662  0.39265  0.49537  0.57807  0.6469   0.70447  0.75265  0.79297  0.82672  0.85497  0.87861]
    >>> qlmFromME = MomentsFromME(alpha, A, 5)
    >>> print(qlmFromME)
    [0.40636152724912805, 0.45459875335830258, 0.76609389063916078, 1.7219814531396285, 4.8383553122056906]
    >>> std, stm = FluidQueue(Q, Rin, Rout, "stDistr", np.arange(0.,1.1,0.1), "stMoms", 5)
    Final Residual Error for G:  6.661338147750939e-15
    >>> print(std)
    [ 0.31678  0.57513  0.67546  0.74313  0.7949   0.83595  0.8688   0.89511  0.91618  0.93304  0.94652]
    >>> print(stm)
    [0.23252260678168307, 0.20068950224745297, 0.26683979953080317, 0.47402477308795288, 1.0522657070139314]
    >>> betap, Bp = FluidQueue(Q, Rin, Rout, "stDistrPH")
    Final Residual Error for G:  6.661338147750939e-15
    >>> print(betap)
    [[ 0.24893  0.0304   0.06859  0.02341  0.21436  0.       0.       0.       0.       0.       0.       0.09753]]
    >>> print(Bp)
    [[ -21.23238    2.489      2.         0.         4.         0.         0.         0.         1.         0.         2.         0.     ]
     [  72.6001  -135.85862    0.         2.         0.         4.         0.         0.         0.         1.         0.         2.     ]
     [   6.         0.       -29.07746    0.82967    5.         0.         3.         0.         7.         0.         4.         0.     ]
     [   0.         6.        24.20003  -67.28621    0.         5.         0.         3.         0.         7.         0.         4.     ]
     [   1.         0.         3.         0.        -4.         0.         0.         0.         0.         0.         0.         0.     ]
     [   0.         1.         0.         3.         0.        -4.         0.         0.         0.         0.         0.         0.     ]
     [   0.         0.         0.         0.         0.         0.        -8.         0.         3.         0.         5.         0.     ]
     [   0.         0.         0.         0.         0.         0.         0.        -8.         0.         3.         0.         5.     ]
     [   7.         0.         3.         0.         0.         0.         2.         0.       -19.11619    1.2445     1.         0.     ]
     [   0.         7.         0.         3.         0.         0.         0.         2.        36.30005  -76.42931    0.         1.     ]
     [   7.         0.         8.         0.         0.         0.         3.         0.         8.         0.       -30.07746    0.82967]
     [   0.         7.         0.         8.         0.         0.         0.         3.         0.         8.        24.20003  -68.28621]]
    >>> beta, B = FluidQueue(Q, Rin, Rout, "stDistrME")
    Final Residual Error for G:  6.661338147750939e-15
    >>> print(beta)
    [[ 0.       0.       0.07883 -0.08752  0.09534  0.       0.      -0.10977  0.10977 -0.10977  0.10977  0.59657]]
    >>> print(B)
    [[ -22.23238   -9.02433   12.67634  -16.41886  -10.24852   30.26497    0.         0.         0.       -20.59361   20.59361    0.     ]
     [  -7.88555  -14.        12.87627   -9.09699   -5.67827   16.76854    0.         0.         0.       -11.41004   11.41004    0.     ]
     [  -7.20983   15.44329  -18.07746  -16.65305    4.82436   14.3251    -4.95592    5.11088    0.       -13.92489   13.92489    0.     ]
     [  -6.49415   13.91033   45.22483  -88.75656    8.12592   15.24683   -9.32892    9.62061   -6.27132   -6.27132   11.28838    1.25426]
     [  -5.17874    9.94388   41.78795  -50.49679  -22.61806   13.99623   -8.56373    3.07457    0.        -6.90832   11.51386    1.15139]
     [   1.15647    9.24387   33.46571  -40.96853  105.26888 -123.84766  -13.15029    9.28013    0.        -6.422      3.211      8.56267]
     [   1.11426    8.90648   32.24425  -39.47323   97.844   -111.61937  -21.67032   14.09775    0.        -6.1876     3.0938     8.25014]
     [   1.08047    5.28041   29.44999  -31.89819   93.1404  -109.16943   47.14189  -53.61589    0.        -6.         0.        11.     ]
     [   1.08047    3.43993   31.60441  -34.29003   95.74595 -109.16943   47.14189  -45.61589   -8.        -8.         2.        11.     ]
     [   1.08047    3.43993   29.44999  -31.89819   95.74595 -109.16943   44.23284  -42.61589    0.       -16.         0.        13.     ]
     [   1.42039   -0.85454   34.47696  -38.27643  102.6941  -109.16943   44.23284  -45.61589    3.         3.11619  -24.62166   18.50548]
     [   1.42039   -0.85454   28.73186  -31.89819  101.82558 -108.23515   37.44508  -38.61589    0.        -2.08935   60.01318  -57.92383]]
    >>> stdFromPH = CdfFromPH(betap, Bp, np.arange(0.,1.1,0.1))
    >>> print(stdFromPH)
    [ 0.31678  0.57513  0.67546  0.74313  0.7949   0.83595  0.8688   0.89511  0.91618  0.93304  0.94652]
    >>> stmFromME = MomentsFromME(beta, B, 5)
    >>> print(stmFromME)
    [0.23252260678168368, 0.20068950224745336, 0.2668397995308035, 0.47402477308795293, 1.0522657070139294]

