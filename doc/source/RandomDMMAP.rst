butools.dmap.RandomDMMAP
========================

.. currentmodule:: butools.dmap

.. np:function:: RandomDMMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`D = RandomDMMAP(order, types, mean, zeroEntries, maxTrials, prec)`
        * - Mathematica:
          - :code:`D = RandomDMMAP[order, types, mean, zeroEntries, maxTrials, prec]`
        * - Python/Numpy:
          - :code:`D = RandomDMMAP(order, types, mean, zeroEntries, maxTrials, prec)`

    Returns a random discrete Markovian arrival process.

    Parameters
    ----------
    order : int
        The size of the DMAP
    mean : double, optional
        The mean inter-arrival times of the DMMAP
    types : int
        The number of different arrival types
    zeroEntries : int, optional
        The number of zero entries in the D0 and D1 matrices
    maxTrials : int, optional
        The maximum number of trials to find a proper DMMAP 
        (that has an irreducible phase process and none of 
        its parameters is all-zero)
    prec : double, optional
        Numerical precision for checking the irreducibility.
        The default value is 1e-14.

    Returns
    -------
    D : list/cell of matrices of shape(M,M), length(types+1)
        The D0...Dtypes matrices of the DMMAP 

    Notes
    -----
    If it fails, try to increase the 'maxTrials' parameter,
    or/and the 'mean' parameter.

    Examples
    ========
    For Matlab:

    >>> Dm = RandomDMMAP(4, 3, 5.62, 10);
    CheckProbMatrix: the matrix has negative element (precision: 1e-12)!
    CheckDMAPRepresentation: D0 isn't a transient probability matrix!
    CheckProbMatrix: the matrix has negative element (precision: 1e-12)!
    CheckDMAPRepresentation: D0 isn't a transient probability matrix!
    >>> disp(Dm{1});
          0.58008      0.05247     0.011422     0.034628
         0.010303      0.78341     0.018769     0.014111
         0.013438     0.016619      0.84305     0.011715
         0.059279            0    0.0073858      0.73922
    >>> disp(Dm{2});
        0.0036364    0.0064271     0.033468     0.029809
         0.019378     0.008356            0     0.038026
        0.0065461     0.015749     0.015421    0.0064125
                0            0     0.045432     0.006937
    >>> disp(Dm{3});
         0.048991     0.045334     0.018459     0.023063
         0.042294      0.00432    0.0039673     0.002284
        0.0010788   0.00023017     0.012396            0
                0     0.040277            0            0
    >>> disp(Dm{4});
         0.026583     0.021537     0.064096            0
          0.01306     0.038122     0.003036   0.00056441
         0.017399     0.014873    0.0088295     0.016249
                0     0.050186      0.04298     0.008303
    >>> m = MarginalMomentsFromDMMAP(Dm, 1);
    >>> disp(m);
             5.62

    For Mathematica:

    >>> Dm = RandomDMMAP[4, 3, 5.62, 10];
    >>> Print[Dm[[1]]];
    {{0.5936500195002852, 0.0071711070249426, 0.017673758024475288, 0.05144659384552179},
     {0.020154515280230965, 0.8411554095030485, 0.005249104682873274, 0.},
     {0.015617402627286933, 0.012537161360729446, 0.8352935249678363, 0.006411658273261797},
     {0.013153759039417926, 0.02868436681898334, 0.028903946905958065, 0.6538234601251783}}
    >>> Print[Dm[[2]]];
    {{0.08519495111289338, 0., 0.018628042747125154, 0.},
     {0.002711391271515091, 0.010446336642803812, 0.01459469095042103, 0.01866365592332515},
     {0., 0.011750889045204358, 0.022017512911504132, 0.012160658337042149},
     {0.03341099177641968, 0.01984209195647928, 0.009966627127729678, 0.04171857977172423}}
    >>> Print[Dm[[3]]];
    {{0.041825271032060074, 0.07761370787172987, 0.08109386986023437, 0.},
     {0.01770089213830925, 0.011074775342567164, 0.017685427551298035, 0.0036482985602606767},
     {0.014365191362523854, 0.015410717929132125, 0.015650791704003995, 0.0195262187884841},
     {0.010576404991469165, 0.03762430513140506, 0.02843564618548617, 0.002495652549692818}}
    >>> Print[Dm[[4]]];
    {{0., 0., 0.025702678980732353, 0.},
     {0., 0.013169854463771804, 0.004281152201241524, 0.019464495488333844},
     {0.001199989751644204, 0.007852396960251728, 0.0014020400199298716, 0.00880384596116504},
     {0.03646020660551431, 0.044521395687431596, 0.010382565327110316, 0.}}
    >>> m = MarginalMomentsFromDMMAP[Dm, 1][[1]];
    >>> Print[m];
    5.62

    For Python/Numpy:

    >>> Dm = RandomDMMAP(4, 3, 5.62, 10)
    >>> print(Dm[0])
    [[ 0.29063  0.03984  0.06259  0.07765]
     [ 0.01556  0.84309  0.01363  0.01053]
     [ 0.       0.00355  0.92541  0.     ]
     [ 0.04962  0.       0.       0.36085]]
    >>> print(Dm[1])
    [[ 0.05324  0.07616  0.00057  0.02539]
     [ 0.01533  0.00616  0.00771  0.01641]
     [ 0.00574  0.00682  0.0007   0.00409]
     [ 0.04219  0.       0.       0.16326]]
    >>> print(Dm[2])
    [[ 0.0679   0.00653  0.06991  0.09437]
     [ 0.00512  0.01717  0.0088   0.00242]
     [ 0.00205  0.00353  0.00687  0.01247]
     [ 0.       0.08672  0.10873  0.16826]]
    >>> print(Dm[3])
    [[ 0.0577   0.03984  0.00719  0.03049]
     [ 0.00606  0.0186   0.00132  0.0121 ]
     [ 0.00942  0.00709  0.00834  0.00392]
     [ 0.       0.       0.       0.02038]]
    >>> m = MarginalMomentsFromDMMAP(Dm, 1)[0]
    >>> print(m)
    5.62

