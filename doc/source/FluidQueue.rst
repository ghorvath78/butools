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
 
    For Python/Numpy:
    
    >>> Q = ml.matrix([[-9, 2, 4, 0, 1, 2], [ 6, -25, 5, 3, 7, 4], [ 1, 3, -4, 0, 0, 0], [ 0, 0, 0, -8, 3, 5], [ 7, 3, 0, 2, -13, 1], [ 7, 8, 0, 3, 8, -26]])
    >>> Rin = Diag(ml.matrix([[4, 2, 1, 0, 0, 3]]))
    >>> Rout = Diag(ml.matrix([[6, 2, 0, 0, 3, 2]]))
    >>> qld,qlm = FluidQueue(Q, Rin, Rout, 'qlDistr', np.arange(0,2.5,0.5), 'qlMoms', 5)
    >>> print(qld)
    [ 0.23662289  0.70447058  0.87860959  0.95013803  0.97951884]
    >>> print(qlm)
    [0.406361527249128, 0.45459875335830258, 0.76609389063916067, 1.7219814531396287, 4.8383553122056897]
    >>> alphap,Ap = FluidQueue(Q, Rin, Rout, 'qlDistrPH')
    >>> print(alphap)
    [[ 0.63124379  0.13213332]]
    >>> print(Ap)
    [[ -2.03872923   0.414833  ]
     [ 12.10001665 -21.14310274]]
    >>> print(CdfFromPH(alphap,Ap,np.arange(0,2.5,0.5)))
    [ 0.23662289  0.70447058  0.87860959  0.95013803  0.97951884]
    >>> print(MomentsFromPH(alphap,Ap,5))
    [0.40636152724912783, 0.45459875335830224, 0.76609389063916011, 1.7219814531396267, 4.8383553122056826]
    >>> std, stm = FluidQueue(Q, Rin, Rout, 'stDistr', np.arange(0,1.2,0.2), 'stMoms', 5)
    >>> print(std)
    [ 0.31678184  0.67545872  0.79489564  0.86879661  0.91618177  0.94652186]
    >>> print(stm)
    [0.23252260678168307, 0.20068950224745297, 0.26683979953080317, 0.47402477308795288, 1.0522657070139314]
    >>> beta, B = FluidQueue(Q, Rin, Rout, 'stDistrME')
    >>> print(beta)
    [[ 0.          0.          0.07883143 -0.08751889  0.09533889  0.          0.  -0.10977188  0.10977188 -0.10977188  0.10977188  0.59656673]]
    >>> print(B)
    [[ -22.23237537   -9.02433446   12.67634128  -16.41886435  -10.24851559    30.26497385    0.            0.            0.          -20.59360623    20.59360623    0.        ]
     [  -7.88555401  -14.           12.8762731    -9.09699459   -5.67826671    16.76853511    0.            0.            0.          -11.41004155    11.41004155    0.        ]
     [  -7.20982894   15.44329182  -18.07745846  -16.65304522    4.82436165    14.32509894   -4.95591772    5.11087606    0.          -13.92488739    13.92488739    0.        ]
     [  -6.49415363   13.91033137   45.224834    -88.75656159    8.125917    15.2468254    -9.32892025    9.62061071   -6.2713248    -6.2713248    11.28838464    1.25426496]
     [  -5.17874382    9.94388398   41.78794639  -50.4967871   -22.61805543    13.99623151   -8.56373206    3.07456666    0.           -6.90831657    11.51386095    1.15138609]
     [   1.15646685    9.24387453   33.46570969  -40.96853204  105.2688808  -123.84766338  -13.15028655    9.28012838    0.           -6.42199885     3.21099943    8.56266514]
     [   1.11425727    8.90648478   32.24425354  -39.47323234   97.84400096  -111.61936992  -21.67031769   14.09775208    0.           -6.187604     3.093802      8.25013867]
     [   1.08047373    5.28041166   29.44999431  -31.89819275   93.14039968  -109.16943403   47.14188668  -53.61588779    0.           -6.            0.    11.        ]
     [   1.08047373    3.43992757   31.6044103   -34.29003187   95.74595486  -109.16943403   47.14188668  -45.61588779   -8.           -8.            2.    11.        ]
     [   1.08047373    3.43992757   29.44999431  -31.89819275   95.74595486  -109.16943403   44.23284467  -42.61588779    0.          -16.            0.    13.        ]
     [   1.42038506   -0.85453532   34.47696494  -38.27643041  102.694102  -109.16943403   44.23284467  -45.61588779    3.            3.11618768   -24.62166422   18.50547654]
     [   1.42038506   -0.85453532   28.73185565  -31.89819275  101.82558361  -108.23514548   37.44507997  -38.61588779    0.           -2.08934647    60.01317816  -57.92383169]]
    >>> print(CdfFromME(beta,B,np.arange(0,1.2,0.2)))
    [ 0.31678184  0.67545872  0.79489564  0.86879661  0.91618177  0.94652186]
    >>> print(MomentsFromME(beta,B,5))
    [0.23252260678168368, 0.20068950224745336, 0.2668397995308035, 0.47402477308795293, 1.0522657070139294]
    
