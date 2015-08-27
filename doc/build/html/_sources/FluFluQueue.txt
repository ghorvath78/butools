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
    For Matlab:

    >>> Qin = [-2., 1., 1.; 2., -5., 3.; 4., 0., -4.];
    >>> vRin = [3.,7.,0.];
    >>> Rin = diag(vRin);
    >>> Qout = [-4., 1., 3.; 6., -8., 2.; 3., 7., -10.];
    >>> vRout = [1.,7.,15.];
    >>> Rout = diag(vRout);
    >>> [fld, flm] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'flDistr', (0.:0.1:1.), 'flMoms', 5);
    Final Residual Error for Psi:    1.0339e-15
    >>> disp(fld);
      Columns 1 through 8
           0.3918      0.47163      0.53819      0.59413       0.6415      0.68193      0.71667      0.74673
      Columns 9 through 11
          0.77292      0.79585      0.81605
    >>> disp(flm);
           0.5357       1.0765       3.4298        14.87       81.162
    >>> [alphap, Ap] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'flDistrPH');
    Final Residual Error for Psi:    1.0339e-15
    >>> disp(alphap);
          0.45573      0.15247
    >>> disp(Ap);
          -2.3405      0.53197
          0.92131      -1.2559
    >>> [alpha, A] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'flDistrME');
    Final Residual Error for Psi:    1.0339e-15
    >>> disp(alpha);
         -0.65561       1.2638
    >>> disp(A);
          -2.1425       1.5194
          0.43807      -1.4538
    >>> fldFromPH = CdfFromPH(alphap, Ap, (0.:0.1:1.));
    >>> disp(fldFromPH);
      Columns 1 through 8
           0.3918      0.47163      0.53819      0.59413       0.6415      0.68193      0.71667      0.74673
      Columns 9 through 11
          0.77292      0.79585      0.81605
    >>> flmFromME = MomentsFromME(alpha, A, 5);
    >>> disp(flmFromME);
           0.5357       1.0765       3.4298        14.87       81.162
    >>> [std, stm] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'stDistr', (0.:0.1:1.), 'stMoms', 5);
    Final Residual Error for Psi:    2.0955e-15
    >>> disp(std);
      Columns 1 through 8
          0.29838      0.51911      0.66679      0.76705      0.83596      0.88381      0.91733      0.94097
      Columns 9 through 11
          0.95774      0.96968      0.97821
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
      Columns 1 through 8
          0.29838      0.51911      0.66679      0.76705      0.83596      0.88381      0.91733      0.94097
      Columns 9 through 11
          0.95774      0.96968      0.97821
    >>> stmFromME = MomentsFromME(beta, B, 5);
    >>> disp(stmFromME);
           0.1948      0.11287      0.10069      0.12158      0.18506
    >>> [fld, flm] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'flDistr', (0.:0.1:1.), 'flMoms', 5);
    Final Residual Error for Psi:    1.0339e-15
    >>> disp(fld);
      Columns 1 through 8
          0.64736      0.68913      0.72467      0.75512       0.7814      0.80423      0.82418      0.84172
      Columns 9 through 11
          0.85721      0.87095      0.88319
    >>> disp(flm);
          0.33265      0.68892       2.2198       9.6621        52.81
    >>> [alphap, Ap] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'flDistrPH');
    Final Residual Error for Psi:    1.0339e-15
    >>> disp(alphap);
          0.24142      0.11122
    >>> disp(Ap);
          -2.3405      0.73252
          0.66907      -1.2559
    >>> [alpha, A] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'flDistrME');
    Final Residual Error for Psi:    1.0339e-15
    >>> disp(alpha);
         -0.24261      0.59524
    >>> disp(A);
          -2.1425       1.5194
          0.43807      -1.4538
    >>> fldFromPH = CdfFromPH(alphap, Ap, (0.:0.1:1.));
    >>> disp(fldFromPH);
      Columns 1 through 8
          0.64736      0.68913      0.72467      0.75512       0.7814      0.80423      0.82418      0.84172
      Columns 9 through 11
          0.85721      0.87095      0.88319
    >>> flmFromME = MomentsFromME(alpha, A, 5);
    >>> disp(flmFromME);
          0.33265      0.68892       2.2198       9.6621        52.81
    >>> [std, stm] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'stDistr', (0.:0.1:1.), 'stMoms', 5);
    Final Residual Error for Psi:    2.0955e-15
    >>> disp(std);
      Columns 1 through 8
          0.57864      0.70628      0.79365      0.85412      0.89636      0.92608      0.94712      0.96209
      Columns 9 through 11
          0.97277      0.98041      0.98589
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
      Columns 1 through 8
          0.57864      0.70628      0.79365      0.85412      0.89636      0.92608      0.94712      0.96209
      Columns 9 through 11
          0.97277      0.98041      0.98589
    >>> stmFromME = MomentsFromME(beta, B, 5);
    >>> disp(stmFromME);
          0.12096     0.071546     0.064592      0.07852      0.11997

    For Mathematica:

    >>> Qin = {{-2., 1., 1.},{2., -5., 3.},{4., 0., -4.}};
    >>> vRin = {3.,7.,0.};
    >>> Rin = DiagonalMatrix[vRin];
    >>> Qout = {{-4., 1., 3.},{6., -8., 2.},{3., 7., -10.}};
    >>> vRout = {1.,7.,15.};
    >>> Rout = DiagonalMatrix[vRout];
    >>> {fld, flm} = FluFluQueue[Qin, Rin, Qout, Rout, False, "flDistr", Range[0.,1.,0.1], "flMoms", 5];
    "Final Residual Error for Psi: "7.008282842946301*^-16
    >>> Print[fld];
    {0.391800873246579, 0.4716300609142923, 0.5381929941518719, 0.59413043637473, 0.641503602438246, 0.6819269051118623, 0.716669775887913, 0.7467348114111322, 0.7729177924650865, 0.7958538190537734, 0.8160528083141847}
    >>> Print[flm];
    {0.5357035627676002, 1.0765385576008892, 3.4298090557057153, 14.869885651621978, 81.16233569132378}
    >>> {alphap, Ap} = FluFluQueue[Qin, Rin, Qout, Rout, False, "flDistrPH"];
    "Final Residual Error for Psi: "7.008282842946301*^-16
    >>> Print[alphap];
    {0.4557324623147923, 0.15246666443862872}
    >>> Print[Ap];
    {{-2.3404575203667424, 0.5319666549172487},
     {0.9213094524423505, -1.2559180746181111}}
    >>> {alpha, A} = FluFluQueue[Qin, Rin, Qout, Rout, False, "flDistrME"];
    "Final Residual Error for Psi: "7.008282842946301*^-16
    >>> Print[alpha];
    {-0.6556057500800456, 1.2638048768334667}
    >>> Print[A];
    {{-2.142529166146859, 1.5193788109637358},
     {0.43806809528730795, -1.4538464288379942}}
    >>> fldFromPH = CdfFromPH[alphap, Ap, Range[0.,1.,0.1]];
    >>> Print[fldFromPH];
    {0.391800873246579, 0.4716300609142924, 0.538192994151872, 0.5941304363747301, 0.641503602438246, 0.6819269051118624, 0.716669775887913, 0.7467348114111322, 0.7729177924650866, 0.7958538190537734, 0.8160528083141847}
    >>> flmFromME = MomentsFromME[alpha, A, 5];
    >>> Print[flmFromME];
    {0.5357035627676006, 1.07653855760089, 3.429809055705716, 14.869885651621974, 81.16233569132378}
    >>> {std, stm} = FluFluQueue[Qin, Rin, Qout, Rout, False, "stDistr", Range[0.,1.,0.1], "stMoms", 5];
    "Final Residual Error for Psi: "1.4988010832439613*^-15
    >>> Print[std];
    {0.29837541247887833, 0.5191050664654031, 0.6667930660737641, 0.7670494936141816, 0.8359575279199422, 0.8838138562939544, 0.9173341131830883, 0.9409745470882221, 0.9577383239502405, 0.9696768240390095, 0.9782073950108067}
    >>> Print[stm];
    {0.19480129555185477, 0.11287108072563842, 0.100686427731073, 0.12157602438536942, 0.18505841424939884}
    >>> {betap, Bp} = FluFluQueue[Qin, Rin, Qout, Rout, False, "stDistrPH"];
    "Final Residual Error for Psi: "1.4988010832439613*^-15
    >>> Print[betap];
    {0.45285051259891584, 0.24877407492220602}
    >>> Print[Bp];
    {{-5.4973439722064095, 0.8367525984737355},
     {1.4491661670964533, -3.7914265223267773}}
    >>> {beta, B} = FluFluQueue[Qin, Rin, Qout, Rout, False, "stDistrME"];
    "Final Residual Error for Psi: "1.4988010832439613*^-15
    >>> Print[beta];
    {0.18114020503956624, 0.5204843824815555}
    >>> Print[B];
    {{-6.366843672464282, 1.6656012008589829},
     {-0.6164326032040903, -2.9219268220689054}}
    >>> stdFromPH = CdfFromPH[betap, Bp, Range[0.,1.,0.1]];
    >>> Print[stdFromPH];
    {0.2983754124788781, 0.519105066465403, 0.666793066073764, 0.7670494936141815, 0.8359575279199422, 0.8838138562939544, 0.9173341131830883, 0.9409745470882221, 0.9577383239502405, 0.9696768240390093, 0.9782073950108067}
    >>> stmFromME = MomentsFromME[beta, B, 5];
    >>> Print[stmFromME];
    {0.19480129555185485, 0.11287108072563846, 0.10068642773107304, 0.12157602438536946, 0.18505841424939895}
    >>> {fld, flm} = FluFluQueue[Qin, Rin, Qout, Rout, True, "flDistr", Range[0.,1.,0.1], "flMoms", 5];
    "Final Residual Error for Psi: "7.008282842946301*^-16
    >>> Print[fld];
    {0.6473632738053204, 0.6891320433434904, 0.7246660281915851, 0.7551190366324814, 0.7814006612465998, 0.8042313747552501, 0.8241848657852636, 0.8417205998183797, 0.8572088895119173, 0.8709502223792169, 0.8831901837025878}
    >>> Print[flm];
    {0.33264746858870425, 0.6889173702248804, 2.219774758007364, 9.662108525950131, 52.80956347845118}
    >>> {alphap, Ap} = FluFluQueue[Qin, Rin, Qout, Rout, True, "flDistrPH"];
    "Final Residual Error for Psi: "7.008282842946301*^-16
    >>> Print[alphap];
    {0.24141936087976823, 0.11121736531491135}
    >>> Print[Ap];
    {{-2.3404575203667424, 0.7325208369425898},
     {0.6690675306998952, -1.2559180746181111}}
    >>> {alpha, A} = FluFluQueue[Qin, Rin, Qout, Rout, True, "flDistrME"];
    "Final Residual Error for Psi: "7.008282842946301*^-16
    >>> Print[alpha];
    {-0.24260693822391072, 0.5952436644185902}
    >>> Print[A];
    {{-2.142529166146859, 1.5193788109637358},
     {0.43806809528730795, -1.4538464288379942}}
    >>> fldFromPH = CdfFromPH[alphap, Ap, Range[0.,1.,0.1]];
    >>> Print[fldFromPH];
    {0.6473632738053204, 0.6891320433434904, 0.7246660281915852, 0.7551190366324814, 0.7814006612465998, 0.8042313747552501, 0.8241848657852637, 0.8417205998183797, 0.8572088895119174, 0.8709502223792169, 0.883190183702588}
    >>> flmFromME = MomentsFromME[alpha, A, 5];
    >>> Print[flmFromME];
    {0.3326474685887042, 0.6889173702248805, 2.219774758007364, 9.662108525950126, 52.80956347845115}
    >>> {std, stm} = FluFluQueue[Qin, Rin, Qout, Rout, True, "stDistr", Range[0.,1.,0.1], "stMoms", 5];
    "Final Residual Error for Psi: "1.4988010832439613*^-15
    >>> Print[std];
    {0.5786381632358244, 0.7062784286059924, 0.7936484851974519, 0.8541200746555163, 0.8963590382366593, 0.9260821681510065, 0.9471221726585012, 0.9620853669984744, 0.9727657494367457, 0.9804107560683182, 0.9858949950823395}
    >>> Print[stm];
    {0.12096271585043802, 0.07154626102115877, 0.0645920699157504, 0.07852033701278222, 0.11996554734440769}
    >>> {betap, Bp} = FluFluQueue[Qin, Rin, Qout, Rout, True, "stDistrPH"];
    "Final Residual Error for Psi: "1.4988010832439613*^-15
    >>> Print[betap];
    {0.2383068326039564, 0.18305500416021936}
    >>> Print[Bp];
    {{-5.4973439722064095, 0.8367525984737355},
     {1.4491661670964533, -3.7914265223267773}}
    >>> {beta, B} = FluFluQueue[Qin, Rin, Qout, Rout, True, "stDistrME"];
    "Final Residual Error for Psi: "1.4988010832439613*^-15
    >>> Print[beta];
    {-0.12203666944014618, 0.5433985062043218}
    >>> Print[B];
    {{-5.18601418644967, 2.4838988377244475},
     {0.6629779802695965, -4.102756308083516}}
    >>> stdFromPH = CdfFromPH[betap, Bp, Range[0.,1.,0.1]];
    >>> Print[stdFromPH];
    {0.5786381632358242, 0.7062784286059923, 0.7936484851974519, 0.8541200746555163, 0.8963590382366593, 0.9260821681510065, 0.9471221726585012, 0.9620853669984744, 0.9727657494367457, 0.9804107560683182, 0.9858949950823395}
    >>> stmFromME = MomentsFromME[beta, B, 5];
    >>> Print[stmFromME];
    {0.12096271585043804, 0.07154626102115878, 0.06459206991575042, 0.07852033701278228, 0.11996554734440777}

    For Python/Numpy:

    >>> Qin = ml.matrix([[-2., 1., 1.],[2., -5., 3.],[4., 0., -4.]])
    >>> vRin = ml.matrix([[3.,7.,0.]])
    >>> Rin = Diag(vRin)
    >>> Qout = ml.matrix([[-4., 1., 3.],[6., -8., 2.],[3., 7., -10.]])
    >>> vRout = ml.matrix([[1.,7.,15.]])
    >>> Rout = Diag(vRout)
    >>> fld, flm = FluFluQueue(Qin, Rin, Qout, Rout, False, "flDistr", np.arange(0.,1.1,0.1), "flMoms", 5)
    Final Residual Error for G:  1.2975731600306517e-15
    >>> print(fld)
    [ 0.3918   0.47163  0.53819  0.59413  0.6415   0.68193  0.71667  0.74673  0.77292  0.79585  0.81605]
    >>> print(flm)
    [0.53570356276760078, 1.0765385576008903, 3.4298090557057193, 14.869885651621992, 81.162335691323889]
    >>> alphap, Ap = FluFluQueue(Qin, Rin, Qout, Rout, False, "flDistrPH")
    Final Residual Error for G:  1.2975731600306517e-15
    >>> print(alphap)
    [[ 0.45573  0.15247]]
    >>> print(Ap)
    [[-2.34046  0.53197]
     [ 0.92131 -1.25592]]
    >>> alpha, A = FluFluQueue(Qin, Rin, Qout, Rout, False, "flDistrME")
    Final Residual Error for G:  1.2975731600306517e-15
    >>> print(alpha)
    [[-0.65561  1.2638 ]]
    >>> print(A)
    [[-2.14253  1.51938]
     [ 0.43807 -1.45385]]
    >>> fldFromPH = CdfFromPH(alphap, Ap, np.arange(0.,1.1,0.1))
    >>> print(fldFromPH)
    [ 0.3918   0.47163  0.53819  0.59413  0.6415   0.68193  0.71667  0.74673  0.77292  0.79585  0.81605]
    >>> flmFromME = MomentsFromME(alpha, A, 5)
    >>> print(flmFromME)
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
    >>> fld, flm = FluFluQueue(Qin, Rin, Qout, Rout, True, "flDistr", np.arange(0.,1.1,0.1), "flMoms", 5)
    Final Residual Error for G:  1.2975731600306517e-15
    >>> print(fld)
    [ 0.64736  0.68913  0.72467  0.75512  0.7814   0.80423  0.82418  0.84172  0.85721  0.87095  0.88319]
    >>> print(flm)
    [0.33264746858870425, 0.68891737022488087, 2.2197747580073659, 9.6621085259501385, 52.809563478451246]
    >>> alphap, Ap = FluFluQueue(Qin, Rin, Qout, Rout, True, "flDistrPH")
    Final Residual Error for G:  1.2975731600306517e-15
    >>> print(alphap)
    [[ 0.24142  0.11122]]
    >>> print(Ap)
    [[-2.34046  0.73252]
     [ 0.66907 -1.25592]]
    >>> alpha, A = FluFluQueue(Qin, Rin, Qout, Rout, True, "flDistrME")
    Final Residual Error for G:  1.2975731600306517e-15
    >>> print(alpha)
    [[-0.24261  0.59524]]
    >>> print(A)
    [[-2.14253  1.51938]
     [ 0.43807 -1.45385]]
    >>> fldFromPH = CdfFromPH(alphap, Ap, np.arange(0.,1.1,0.1))
    >>> print(fldFromPH)
    [ 0.64736  0.68913  0.72467  0.75512  0.7814   0.80423  0.82418  0.84172  0.85721  0.87095  0.88319]
    >>> flmFromME = MomentsFromME(alpha, A, 5)
    >>> print(flmFromME)
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

