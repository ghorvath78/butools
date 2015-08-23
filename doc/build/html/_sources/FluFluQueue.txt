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
    For Matlab:

    >>> Qin = [-2., 1., 1.; 2., -5., 3.; 4., 0., -4.];
    >>> vRin = [3.,7.,0.];
    >>> Rin = diag(vRin);
    >>> Qout = [-4., 1., 3.; 6., -8., 2.; 3., 7., -10.];
    >>> vRout = [1.,7.,15.];
    >>> Rout = diag(vRout);
    >>> [qld, qlm] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'qlDistr', (0.:0.1:1.), 'qlMoms', 5);
    Final Residual Error for Psi:    1.0339e-15
    >>> disp(qld);
      Columns 1 through 6
           0.3918      0.47163      0.53819      0.59413       0.6415      0.68193
      Columns 7 through 11
          0.71667      0.74673      0.77292      0.79585      0.81605
    >>> disp(qlm);
           0.5357       1.0765       3.4298        14.87       81.162
    >>> [alphap, Ap] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'qlDistrPH');
    Final Residual Error for Psi:    1.0339e-15
    >>> disp(alphap);
          0.45573      0.15247
    >>> disp(Ap);
          -2.3405      0.53197
          0.92131      -1.2559
    >>> [alpha, A] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'qlDistrME');
    Final Residual Error for Psi:    1.0339e-15
    >>> disp(alpha);
         -0.65561       1.2638
    >>> disp(A);
          -2.1425       1.5194
          0.43807      -1.4538
    >>> qldFromPH = CdfFromPH(alphap, Ap, (0.:0.1:1.));
    >>> disp(qldFromPH);
      Columns 1 through 6
           0.3918      0.47163      0.53819      0.59413       0.6415      0.68193
      Columns 7 through 11
          0.71667      0.74673      0.77292      0.79585      0.81605
    >>> qlmFromME = MomentsFromME(alpha, A, 5);
    >>> disp(qlmFromME);
           0.5357       1.0765       3.4298        14.87       81.162
    >>> [std, stm] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'stDistr', (0.:0.1:1.), 'stMoms', 5);
    Final Residual Error for Psi:    2.0955e-15
    >>> disp(std);
      Columns 1 through 6
          0.29838      0.51911      0.66679      0.76705      0.83596      0.88381
      Columns 7 through 11
          0.91733      0.94097      0.95774      0.96968      0.97821
    >>> disp(stm);
           0.1948      0.11287      0.10069      0.12158      0.18506
    >>> [betap, Bp] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'stDistrPH');
    Final Residual Error for Psi:    2.0955e-15
    >>> disp(betap);
          0.45285      0.24877
    >>> disp(Bp);
          -5.4973      0.83675
           1.4492      -3.7914
    >>> [beta, B] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'stDistrME');
    Final Residual Error for Psi:    2.0955e-15
    >>> disp(beta);
          0.18114      0.52048
    >>> disp(B);
          -6.3668       1.6656
         -0.61643      -2.9219
    >>> stdFromPH = CdfFromPH(betap, Bp, (0.:0.1:1.));
    >>> disp(stdFromPH);
      Columns 1 through 6
          0.29838      0.51911      0.66679      0.76705      0.83596      0.88381
      Columns 7 through 11
          0.91733      0.94097      0.95774      0.96968      0.97821
    >>> stmFromME = MomentsFromME(beta, B, 5);
    >>> disp(stmFromME);
           0.1948      0.11287      0.10069      0.12158      0.18506
    >>> [qld, qlm] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'qlDistr', (0.:0.1:1.), 'qlMoms', 5);
    Final Residual Error for Psi:    1.0339e-15
    >>> disp(qld);
      Columns 1 through 6
          0.64736      0.68913      0.72467      0.75512       0.7814      0.80423
      Columns 7 through 11
          0.82418      0.84172      0.85721      0.87095      0.88319
    >>> disp(qlm);
          0.33265      0.68892       2.2198       9.6621        52.81
    >>> [alphap, Ap] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'qlDistrPH');
    Final Residual Error for Psi:    1.0339e-15
    >>> disp(alphap);
          0.24142      0.11122
    >>> disp(Ap);
          -2.3405      0.73252
          0.66907      -1.2559
    >>> [alpha, A] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'qlDistrME');
    Final Residual Error for Psi:    1.0339e-15
    >>> disp(alpha);
         -0.24261      0.59524
    >>> disp(A);
          -2.1425       1.5194
          0.43807      -1.4538
    >>> qldFromPH = CdfFromPH(alphap, Ap, (0.:0.1:1.));
    >>> disp(qldFromPH);
      Columns 1 through 6
          0.64736      0.68913      0.72467      0.75512       0.7814      0.80423
      Columns 7 through 11
          0.82418      0.84172      0.85721      0.87095      0.88319
    >>> qlmFromME = MomentsFromME(alpha, A, 5);
    >>> disp(qlmFromME);
          0.33265      0.68892       2.2198       9.6621        52.81
    >>> [std, stm] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'stDistr', (0.:0.1:1.), 'stMoms', 5);
    Final Residual Error for Psi:    2.0955e-15
    >>> disp(std);
      Columns 1 through 6
          0.57864      0.70628      0.79365      0.85412      0.89636      0.92608
      Columns 7 through 11
          0.94712      0.96209      0.97277      0.98041      0.98589
    >>> disp(stm);
          0.12096     0.071546     0.064592      0.07852      0.11997
    >>> [betap, Bp] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'stDistrPH');
    Final Residual Error for Psi:    2.0955e-15
    >>> disp(betap);
          0.23831      0.18306
    >>> disp(Bp);
          -5.4973      0.83675
           1.4492      -3.7914
    >>> [beta, B] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'stDistrME');
    Final Residual Error for Psi:    2.0955e-15
    >>> disp(beta);
         -0.12204       0.5434
    >>> disp(B);
           -5.186       2.4839
          0.66298      -4.1028
    >>> stdFromPH = CdfFromPH(betap, Bp, (0.:0.1:1.));
    >>> disp(stdFromPH);
      Columns 1 through 6
          0.57864      0.70628      0.79365      0.85412      0.89636      0.92608
      Columns 7 through 11
          0.94712      0.96209      0.97277      0.98041      0.98589
    >>> stmFromME = MomentsFromME(beta, B, 5);
    >>> disp(stmFromME);
          0.12096     0.071546     0.064592      0.07852      0.11997

    For Mathematica:

    
    For Python/Numpy:

    >>> Qin = ml.matrix([[-2., 1., 1.],[2., -5., 3.],[4., 0., -4.]])
    >>> vRin = ml.matrix([[3.,7.,0.]])
    >>> Rin = Diag(vRin)
    >>> Qout = ml.matrix([[-4., 1., 3.],[6., -8., 2.],[3., 7., -10.]])
    >>> vRout = ml.matrix([[1.,7.,15.]])
    >>> Rout = Diag(vRout)
    >>> qld, qlm = FluFluQueue(Qin, Rin, Qout, Rout, False, "qlDistr", np.arange(0.,1.1,0.1), "qlMoms", 5)
    Final Residual Error for G:  1.2975731600306517e-15
    >>> print(qld)
    [ 0.3918   0.47163  0.53819  0.59413  0.6415   0.68193  0.71667  0.74673  0.77292  0.79585  0.81605]
    >>> print(qlm)
    [0.53570356276760078, 1.0765385576008903, 3.4298090557057193, 14.869885651621992, 81.162335691323889]
    >>> alphap, Ap = FluFluQueue(Qin, Rin, Qout, Rout, False, "qlDistrPH")
    Final Residual Error for G:  1.2975731600306517e-15
    >>> print(alphap)
    [[ 0.45573  0.15247]]
    >>> print(Ap)
    [[-2.34046  0.53197]
     [ 0.92131 -1.25592]]
    >>> alpha, A = FluFluQueue(Qin, Rin, Qout, Rout, False, "qlDistrME")
    Final Residual Error for G:  1.2975731600306517e-15
    >>> print(alpha)
    [[-0.65561  1.2638 ]]
    >>> print(A)
    [[-2.14253  1.51938]
     [ 0.43807 -1.45385]]
    >>> qldFromPH = CdfFromPH(alphap, Ap, np.arange(0.,1.1,0.1))
    >>> print(qldFromPH)
    [ 0.3918   0.47163  0.53819  0.59413  0.6415   0.68193  0.71667  0.74673  0.77292  0.79585  0.81605]
    >>> qlmFromME = MomentsFromME(alpha, A, 5)
    >>> print(qlmFromME)
    [0.53570356276760067, 1.0765385576008901, 3.4298090557057166, 14.869885651621992, 81.162335691323918]
    >>> std, stm = FluFluQueue(Qin, Rin, Qout, Rout, False, "stDistr", np.arange(0.,1.1,0.1), "stMoms", 5)
    Final Residual Error for G:  1.8457457784393227e-15
    >>> print(std)
    [[ 0.29838  0.51911  0.66679  0.76705  0.83596  0.88381  0.91733  0.94097  0.95774  0.96968  0.97821]]
    >>> print(stm)
    [0.19480129555185471, 0.1128710807256384, 0.10068642773107296, 0.12157602438536934, 0.18505841424939878]
    >>> betap, Bp = FluFluQueue(Qin, Rin, Qout, Rout, False, "stDistrPH")
    Final Residual Error for G:  1.8457457784393227e-15
    >>> print(betap)
    [[ 0.45285  0.24877]]
    >>> print(Bp)
    [[-5.49734  0.83675]
     [ 1.44917 -3.79143]]
    >>> beta, B = FluFluQueue(Qin, Rin, Qout, Rout, False, "stDistrME")
    Final Residual Error for G:  1.8457457784393227e-15
    >>> print(beta)
    [[ 0.18114  0.52048]]
    >>> print(B)
    [[-6.36684  1.6656 ]
     [-0.61643 -2.92193]]
    >>> stdFromPH = CdfFromPH(betap, Bp, np.arange(0.,1.1,0.1))
    >>> print(stdFromPH)
    [ 0.29838  0.51911  0.66679  0.76705  0.83596  0.88381  0.91733  0.94097  0.95774  0.96968  0.97821]
    >>> stmFromME = MomentsFromME(beta, B, 5)
    >>> print(stmFromME)
    [0.19480129555185474, 0.11287108072563844, 0.10068642773107304, 0.12157602438536946, 0.18505841424939901]
    >>> qld, qlm = FluFluQueue(Qin, Rin, Qout, Rout, True, "qlDistr", np.arange(0.,1.1,0.1), "qlMoms", 5)
    Final Residual Error for G:  1.2975731600306517e-15
    >>> print(qld)
    [ 0.64736  0.68913  0.72467  0.75512  0.7814   0.80423  0.82418  0.84172  0.85721  0.87095  0.88319]
    >>> print(qlm)
    [0.33264746858870425, 0.68891737022488087, 2.2197747580073659, 9.6621085259501385, 52.809563478451246]
    >>> alphap, Ap = FluFluQueue(Qin, Rin, Qout, Rout, True, "qlDistrPH")
    Final Residual Error for G:  1.2975731600306517e-15
    >>> print(alphap)
    [[ 0.24142  0.11122]]
    >>> print(Ap)
    [[-2.34046  0.73252]
     [ 0.66907 -1.25592]]
    >>> alpha, A = FluFluQueue(Qin, Rin, Qout, Rout, True, "qlDistrME")
    Final Residual Error for G:  1.2975731600306517e-15
    >>> print(alpha)
    [[-0.24261  0.59524]]
    >>> print(A)
    [[-2.14253  1.51938]
     [ 0.43807 -1.45385]]
    >>> qldFromPH = CdfFromPH(alphap, Ap, np.arange(0.,1.1,0.1))
    >>> print(qldFromPH)
    [ 0.64736  0.68913  0.72467  0.75512  0.7814   0.80423  0.82418  0.84172  0.85721  0.87095  0.88319]
    >>> qlmFromME = MomentsFromME(alpha, A, 5)
    >>> print(qlmFromME)
    [0.3326474685887042, 0.68891737022488064, 2.219774758007365, 9.6621085259501367, 52.809563478451238]
    >>> std, stm = FluFluQueue(Qin, Rin, Qout, Rout, True, "stDistr", np.arange(0.,1.1,0.1), "stMoms", 5)
    Final Residual Error for G:  1.8457457784393227e-15
    >>> print(std)
    [[ 0.57864  0.70628  0.79365  0.85412  0.89636  0.92608  0.94712  0.96209  0.97277  0.98041  0.98589]]
    >>> print(stm)
    [0.12096271585043798, 0.071546261021158727, 0.064592069915750366, 0.078520337012782174, 0.11996554734440762]
    >>> betap, Bp = FluFluQueue(Qin, Rin, Qout, Rout, True, "stDistrPH")
    Final Residual Error for G:  1.8457457784393227e-15
    >>> print(betap)
    [[ 0.23831  0.18306]]
    >>> print(Bp)
    [[-5.49734  0.83675]
     [ 1.44917 -3.79143]]
    >>> beta, B = FluFluQueue(Qin, Rin, Qout, Rout, True, "stDistrME")
    Final Residual Error for G:  1.8457457784393227e-15
    >>> print(beta)
    [[-0.12204  0.5434 ]]
    >>> print(B)
    [[-5.18601  2.4839 ]
     [ 0.66298 -4.10276]]
    >>> stdFromPH = CdfFromPH(betap, Bp, np.arange(0.,1.1,0.1))
    >>> print(stdFromPH)
    [ 0.57864  0.70628  0.79365  0.85412  0.89636  0.92608  0.94712  0.96209  0.97277  0.98041  0.98589]
    >>> stmFromME = MomentsFromME(beta, B, 5)
    >>> print(stmFromME)
    [0.12096271585043798, 0.071546261021158755, 0.064592069915750394, 0.078520337012782215, 0.11996554734440766]

