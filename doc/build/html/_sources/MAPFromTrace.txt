butools.fitting.MAPFromTrace
============================

.. currentmodule:: butools.fitting

.. np:function:: MAPFromTrace

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[D0, D1, logli] = MAPFromTrace(trace, orders, maxIter, stopCond, initial, result)`
        * - Mathematica:
          - :code:`{D0, D1, logli} = MAPFromTrace[trace, orders, maxIter, stopCond, initial, result]`
        * - Python/Numpy:
          - :code:`D0, D1, logli = MAPFromTrace(trace, orders, maxIter, stopCond, initial, result)`

    Performs MAP fitting using the EM algorithm (ErCHMM, 
    [1]_, [2]_).
    
    Parameters
    ----------
    trace : column vector, length K
        The samples of the trace
    orders : list of int, length(N), or int
        The length of the list determines the number of 
        Erlang branches to use in the fitting method.
        The entries of the list are the orders of the 
        Erlang distributions. If this parameter is a 
        single integer, all possible branch number - order
        combinations are tested where the total number of 
        states is "orders".
    maxIter : int, optional
        Maximum number of iterations. The default value is
        200
    stopCond : double, optional
        The algorithm stops if the relative improvement of
        the log likelihood falls below stopCond. The 
        default value is 1e-7
    initial : tuple of a vector and a matrix, shape(N,N), optional
        The rate parameters of the Erlang distributions 
        and the branch transition probability matrix to be
        used initially. If not given, a default initial 
        guess is determined and the algorithm starts from 
        there.
    result : {"vecmat", "matmat"}, optional
        The result can be returned two ways. If "matmat" is
        selected, the result is returned in the classical
        representation of MAPs, thus the D0 and D1 matrices.
        If "vecmat" is selected, the rate parameters of the
        Erlang branches and the branch transition probability
        matrix are returned. The default value is "matmat"

    Returns
    -------
    (D0, D1) : tuple of matrix, shape (M,M) and matrix, shape (M,M)
        If the "matmat" result format is chosen, the function
        returns the D0 and D1 matrices of the MAP
    (lambda, P) : tuple of vector, length N and matrix, shape (M,M)
        If the "vecmat" result format is chosen, the function
        returns the vector of the Erlang rate parameters of 
        the branches and the branch transition probability 
        matrix
    logli : double
        The log-likelihood divided by the trace length
        
    Notes
    -----
    This procedure is quite slow in the supported 
    mathematical frameworks. If the maximum speed is
    needed, please use the multi-core optimized c++
    implementation called SPEM-FIT_.

    .. _SPEM-FIT: https://bitbucket.org/ghorvath78/spemfit

    References
    ----------
    .. [1] Okamura, Hiroyuki, and Tadashi Dohi. Faster 
           maximum likelihood estimation algorithms for 
           Markovian arrival processes. Quantitative 
           Evaluation of Systems, 2009. QEST'09. Sixth 
           International Conference on the. IEEE, 2009.
    
    .. [2] Horváth, Gábor, and Hiroyuki Okamura. A Fast EM
           Algorithm for Fitting Marked Markovian Arrival 
           Processes with a New Special Structure. Computer
           Performance Engineering. Springer Berlin 
           Heidelberg, 2013. 119-133.
    
    Examples
    ========
    For Matlab:

    >>> tr = dlmread('/home/gabor/github/butools/test/data/bctrace.iat');
    >>> tr = tr(1:10000);
    >>> [D0, D1] = MAPFromTrace(tr, 5);
    Trying orders 1,4...
    Num of iterations: 8, logli: 4.99752
    Num of iterations: 15, logli: 4.99807
    Num of iterations: 22, logli: 4.9981
    Num of iterations: 26, logli: 4.99811
    EM algorithm terminated. (orders=1,4)
    Trying orders 2,3...
    Num of iterations: 8, logli: 4.93956
    Num of iterations: 14, logli: 4.95626
    Num of iterations: 21, logli: 4.95641
    Num of iterations: 21, logli: 4.95641
    EM algorithm terminated. (orders=2,3)
    Trying orders 1,1,3...
    Num of iterations: 8, logli: 5.04759
    Num of iterations: 14, logli: 5.09151
    Num of iterations: 21, logli: 5.09701
    Num of iterations: 27, logli: 5.09814
    Num of iterations: 34, logli: 5.09943
    Num of iterations: 41, logli: 5.10101
    Num of iterations: 48, logli: 5.10369
    Num of iterations: 55, logli: 5.11054
    Num of iterations: 62, logli: 5.12063
    Num of iterations: 69, logli: 5.12327
    Num of iterations: 76, logli: 5.12353
    Num of iterations: 83, logli: 5.12366
    Num of iterations: 90, logli: 5.12379
    Num of iterations: 97, logli: 5.1239
    Num of iterations: 104, logli: 5.12394
    Num of iterations: 111, logli: 5.12396
    Num of iterations: 115, logli: 5.12396
    EM algorithm terminated. (orders=1,1,3)
    Trying orders 1,2,2...
    Num of iterations: 7, logli: 5.01977
    Num of iterations: 13, logli: 5.04362
    Num of iterations: 19, logli: 5.07693
    Num of iterations: 26, logli: 5.11112
    Num of iterations: 32, logli: 5.11443
    Num of iterations: 39, logli: 5.11519
    Num of iterations: 46, logli: 5.11543
    Num of iterations: 53, logli: 5.11548
    Num of iterations: 59, logli: 5.11549
    Num of iterations: 60, logli: 5.11549
    EM algorithm terminated. (orders=1,2,2)
    Trying orders 1,1,1,2...
    Num of iterations: 7, logli: 5.04003
    Num of iterations: 13, logli: 5.07886
    Num of iterations: 19, logli: 5.08357
    Num of iterations: 25, logli: 5.08526
    Num of iterations: 31, logli: 5.08631
    Num of iterations: 37, logli: 5.08695
    Num of iterations: 43, logli: 5.08738
    Num of iterations: 49, logli: 5.08772
    Num of iterations: 55, logli: 5.08807
    Num of iterations: 61, logli: 5.08853
    Num of iterations: 67, logli: 5.08926
    Num of iterations: 73, logli: 5.09065
    Num of iterations: 79, logli: 5.09386
    Num of iterations: 85, logli: 5.10165
    Num of iterations: 91, logli: 5.11062
    Num of iterations: 97, logli: 5.11257
    Num of iterations: 103, logli: 5.11268
    Num of iterations: 107, logli: 5.11268
    EM algorithm terminated. (orders=1,1,1,2)
    Trying orders 1,1,1,1,1...
    Num of iterations: 7, logli: 5.016
    Num of iterations: 13, logli: 5.04173
    Num of iterations: 19, logli: 5.04393
    Num of iterations: 23, logli: 5.04394
    EM algorithm terminated. (orders=1,1,1,1,1)
    Best solution: logli=5.12396, orders=1,1,3
    >>> disp(D0);
          -83.429            0            0            0            0
                0      -718.68            0            0            0
                0            0      -1026.2       1026.2            0
                0            0            0      -1026.2       1026.2
                0            0            0            0      -1026.2
    >>> disp(D1);
           54.149       4.9019       24.379            0            0
           3.3915       665.85       49.439            0            0
                0            0            0            0            0
                0            0            0            0            0
           42.647       96.944       886.57            0            0
    >>> logli = LikelihoodFromTrace(tr, D0, D1);
    >>> disp(logli);
            5.124
    >>> trAcf = LagCorrelationsFromTrace(tr, 10);
    >>> disp(trAcf);
      Columns 1 through 6
          0.18412      0.18159      0.17544      0.19965     0.083228      0.08634
      Columns 7 through 10
            0.095     0.062854      0.06232     0.065922
    >>> mapAcf = LagCorrelationsFromMAP(D0, D1, 10);
    >>> disp(mapAcf);
      Columns 1 through 6
          0.24889      0.17665      0.12882     0.096383     0.073802      0.05765
      Columns 7 through 10
         0.045782      0.03684     0.029952     0.024544
    >>> sqAcf = SquaredDifference(mapAcf, trAcf);
    >>> disp(sqAcf);
         0.023828
    >>> reAcf = RelativeEntropy(mapAcf, trAcf);
    >>> disp(reAcf);
          0.32132

    For Mathematica:

    >>> tr = Flatten[Import["/home/gabor/github/butools/test/data/bctrace.iat","CSV"]];
    >>> tr = tr[[1;;10000]];
    >>> {D0, D1, logli} = MAPFromTrace[tr, 5];
    "Trying orders "{1, 4}"..."
    "Num of iterations: "4", logli: "4.993259505655185
    "Num of iterations: "7", logli: "4.997168671547701
    "Num of iterations: "10", logli: "4.997862551063841
    "Num of iterations: "13", logli: "4.998033171679149
    "Num of iterations: "16", logli: "4.998082939289365
    "Num of iterations: "19", logli: "4.99809905455834
    "Num of iterations: "22", logli: "4.998104743972023
    "Num of iterations: "25", logli: "4.998107006878427
    "Num of iterations: "26", logli: "4.99810745011236
    "EM algorithm terminated. (orders="{1, 4}")"
    "Trying orders "{2, 3}"..."
    "Num of iterations: "4", logli: "4.8663721991321305
    "Num of iterations: "7", logli: "4.927252921982753
    "Num of iterations: "10", logli: "4.951919574051944
    "Num of iterations: "13", logli: "4.956040780340811
    "Num of iterations: "16", logli: "4.9563853575360515
    "Num of iterations: "19", logli: "4.956407398543356
    "Num of iterations: "21", logli: "4.956408840810655
    "EM algorithm terminated. (orders="{2, 3}")"
    "Trying orders "{1, 1, 3}"..."
    "Num of iterations: "4", logli: "5.015729326700379
    "Num of iterations: "7", logli: "5.039376105933879
    "Num of iterations: "10", logli: "5.065115109858913
    "Num of iterations: "13", logli: "5.087433305519177
    "Num of iterations: "16", logli: "5.095165218779676
    "Num of iterations: "19", logli: "5.096542942667858
    "Num of iterations: "22", logli: "5.097217245266611
    "Num of iterations: "25", logli: "5.097783794710717
    "Num of iterations: "28", logli: "5.098316830378618
    "Num of iterations: "31", logli: "5.098855792506172
    "Num of iterations: "34", logli: "5.099426493035749
    "Num of iterations: "37", logli: "5.100047995275797
    "Num of iterations: "40", logli: "5.100746872797194
    "Num of iterations: "43", logli: "5.101583155259961
    "Num of iterations: "46", logli: "5.1026877401922155
    "Num of iterations: "49", logli: "5.104304334981796
    "Num of iterations: "52", logli: "5.106804300310496
    "Num of iterations: "55", logli: "5.110538394539443
    "Num of iterations: "58", logli: "5.115263627401993
    "Num of iterations: "61", logli: "5.119570638369141
    "Num of iterations: "64", logli: "5.122077926685482
    "Num of iterations: "67", logli: "5.1230372177203165
    "Num of iterations: "70", logli: "5.123338507017323
    "Num of iterations: "73", logli: "5.123456224814315
    "Num of iterations: "76", logli: "5.123528798414036
    "Num of iterations: "79", logli: "5.123585724468395
    "Num of iterations: "82", logli: "5.123637966952736
    "Num of iterations: "85", logli: "5.12369334181503
    "Num of iterations: "88", logli: "5.123752973774425
    "Num of iterations: "91", logli: "5.123811099333816
    "Num of iterations: "94", logli: "5.1238604907978615
    "Num of iterations: "97", logli: "5.12389743186796
    "Num of iterations: "100", logli: "5.123922405926877
    "Num of iterations: "103", logli: "5.123938104622403
    "Num of iterations: "106", logli: "5.123947499602905
    "Num of iterations: "109", logli: "5.123952946957953
    "Num of iterations: "112", logli: "5.12395604420125
    "Num of iterations: "115", logli: "5.123957785028128
    "Num of iterations: "115", logli: "5.123957785028128
    "EM algorithm terminated. (orders="{1, 1, 3}")"
    "Trying orders "{1, 2, 2}"..."
    "Num of iterations: "4", logli: "5.0122947542928555
    "Num of iterations: "7", logli: "5.019767395602894
    "Num of iterations: "10", logli: "5.02976587459673
    "Num of iterations: "13", logli: "5.043619440764079
    "Num of iterations: "16", logli: "5.059183637783567
    "Num of iterations: "19", logli: "5.07692712660751
    "Num of iterations: "22", logli: "5.0973408269250005
    "Num of iterations: "25", logli: "5.10934742598377
    "Num of iterations: "28", logli: "5.113031789500797
    "Num of iterations: "31", logli: "5.114217422911263
    "Num of iterations: "34", logli: "5.1147363762308915
    "Num of iterations: "37", logli: "5.115045193666414
    "Num of iterations: "40", logli: "5.115247596938485
    "Num of iterations: "43", logli: "5.115368318400923
    "Num of iterations: "46", logli: "5.115433152969141
    "Num of iterations: "49", logli: "5.115465572463486
    "Num of iterations: "52", logli: "5.1154810563434925
    "Num of iterations: "55", logli: "5.115488231281802
    "Num of iterations: "58", logli: "5.11549148795029
    "Num of iterations: "60", logli: "5.115492582889059
    "EM algorithm terminated. (orders="{1, 2, 2}")"
    "Trying orders "{1, 1, 1, 2}"..."
    "Num of iterations: "4", logli: "5.014829338214883
    "Num of iterations: "7", logli: "5.040028946911346
    "Num of iterations: "10", logli: "5.064171862323919
    "Num of iterations: "13", logli: "5.078863209647303
    "Num of iterations: "16", logli: "5.082324376375294
    "Num of iterations: "19", logli: "5.083569979450504
    "Num of iterations: "22", logli: "5.084513264289563
    "Num of iterations: "25", logli: "5.08526460887892
    "Num of iterations: "28", logli: "5.0858535498804205
    "Num of iterations: "31", logli: "5.086311241542606
    "Num of iterations: "34", logli: "5.086668013733888
    "Num of iterations: "37", logli: "5.08695043549303
    "Num of iterations: "40", logli: "5.087180397015339
    "Num of iterations: "43", logli: "5.08737551075806
    "Num of iterations: "46", logli: "5.0875500846222375
    "Num of iterations: "49", logli: "5.087716242054033
    "Num of iterations: "52", logli: "5.0878850486431135
    "Num of iterations: "55", logli: "5.088067682179787
    "Num of iterations: "58", logli: "5.088276806171391
    "Num of iterations: "61", logli: "5.088528454903264
    "Num of iterations: "64", logli: "5.088845013008133
    "Num of iterations: "67", logli: "5.089260445776136
    "Num of iterations: "70", logli: "5.089830124804416
    "Num of iterations: "73", logli: "5.090649849553919
    "Num of iterations: "76", logli: "5.091891633656121
    "Num of iterations: "79", logli: "5.093859592976077
    "Num of iterations: "82", logli: "5.097011166206567
    "Num of iterations: "85", logli: "5.1016492083255836
    "Num of iterations: "88", logli: "5.106896654028055
    "Num of iterations: "91", logli: "5.110618528293328
    "Num of iterations: "94", logli: "5.112150973773346
    "Num of iterations: "97", logli: "5.112567107521197
    "Num of iterations: "100", logli: "5.11265838511872
    "Num of iterations: "103", logli: "5.11267690302395
    "Num of iterations: "106", logli: "5.1126806540406635
    "Num of iterations: "107", logli: "5.112681069700427
    "EM algorithm terminated. (orders="{1, 1, 1, 2}")"
    "Trying orders "{1, 1, 1, 1, 1}"..."
    "Num of iterations: "4", logli: "5.001018278123578
    "Num of iterations: "7", logli: "5.016002770232654
    "Num of iterations: "10", logli: "5.031545952575503
    "Num of iterations: "13", logli: "5.041731506932355
    "Num of iterations: "16", logli: "5.043786339264612
    "Num of iterations: "19", logli: "5.04393273379092
    "Num of iterations: "22", logli: "5.043940070674873
    "Num of iterations: "23", logli: "5.043940317578035
    "EM algorithm terminated. (orders="{1, 1, 1, 1, 1}")"
    "Best solution: logli="5.123957785028128", orders="{1, 1, 3}
    >>> Print[D0];
    {{-83.42943388846042, 0, 0, 0, 0},
     {0, -718.6779892150017, 0, 0, 0},
     {0, 0, -1026.160629451268, 1026.160629451268, 0.},
     {0, 0, 0., -1026.160629451268, 1026.160629451268},
     {0, 0, 0., 0., -1026.160629451268}}
    >>> Print[D1];
    {{54.14857308767231, 4.901864224679702, 24.378996576108406, 0., 0.},
     {3.391518630791651, 665.8473464602738, 49.439124123936224, 0., 0.},
     {0., 0., 0., 0., 0.},
     {0., 0., 0., 0., 0.},
     {42.64730201560431, 96.94396093846602, 886.5693664971978, 0., 0.}}
    >>> Print[logli];
    5.123957785028128
    >>> logli = LikelihoodFromTrace[tr, D0, D1];
    >>> Print[logli];
    5.123958173279034
    >>> trAcf = LagCorrelationsFromTrace[tr, 10];
    >>> Print[trAcf];
    {0.18411691802386188, 0.18158522314048947, 0.17543727813332094, 0.19964686019458905, 0.08322774966768201, 0.08633973760499067, 0.0950004809602402, 0.0628536515187086, 0.06232004520613776, 0.06592211464475752}
    >>> mapAcf = LagCorrelationsFromMAP[D0, D1, 10];
    >>> Print[mapAcf];
    {0.248892052255008, 0.17665283034500237, 0.12882046679167725, 0.0963829330310287, 0.07380230626995823, 0.057649876260931265, 0.045781783730242616, 0.03684041958390035, 0.029952264846717966, 0.024544285061666418}
    >>> sqAcf = SquaredDifference[mapAcf, trAcf];
    >>> Print[sqAcf];
    0.023827625772555785
    >>> reAcf = RelativeEntropy[mapAcf, trAcf];
    >>> Print[reAcf];
    0.3213200829124291

    For Python/Numpy:

    >>> tr = np.loadtxt("/home/gabor/github/butools/test/data/bctrace.iat")
    >>> tr = tr[0:10000]
    >>> D0, D1 = MAPFromTrace(tr, 5)
    Trying orders  [1, 4]
    iteration:  10 , logli:  4.99786255106
    iteration:  20 , logli:  4.99810160119
    Num of iterations:  26 , logli:  4.99810745011
    EM algorithm terminated. [1, 4]
    Trying orders  [2, 3]
    iteration:  10 , logli:  4.95191957405
    iteration:  20 , logli:  4.95640838807
    Num of iterations:  21 , logli:  4.95640884081
    EM algorithm terminated. [2, 3]
    Trying orders  [1, 1, 3]
    iteration:  10 , logli:  5.06511510986
    iteration:  20 , logli:  5.09678944434
    iteration:  30 , logli:  5.0986736619
    iteration:  40 , logli:  5.1007468728
    iteration:  50 , logli:  5.10501858379
    iteration:  60 , logli:  5.11829491111
    iteration:  70 , logli:  5.12333850702
    iteration:  80 , logli:  5.12360318453
    iteration:  90 , logli:  5.12379237877
    iteration:  100 , logli:  5.12392240593
    iteration:  110 , logli:  5.123954181
    Num of iterations:  115 , logli:  5.12395778503
    EM algorithm terminated. [1, 1, 3]
    Trying orders  [1, 2, 2]
    iteration:  10 , logli:  5.0297658746
    iteration:  20 , logli:  5.0838748508
    iteration:  30 , logli:  5.11393811996
    iteration:  40 , logli:  5.11524759694
    iteration:  50 , logli:  5.11547207462
    iteration:  60 , logli:  5.11549258289
    Num of iterations:  60 , logli:  5.11549258289
    EM algorithm terminated. [1, 2, 2]
    Trying orders  [1, 1, 1, 2]
    iteration:  10 , logli:  5.06417186232
    iteration:  20 , logli:  5.08390784667
    iteration:  30 , logli:  5.08617133009
    iteration:  40 , logli:  5.08718039702
    iteration:  50 , logli:  5.08777167875
    iteration:  60 , logli:  5.08843863999
    iteration:  70 , logli:  5.0898301248
    iteration:  80 , logli:  5.09475385589
    iteration:  90 , logli:  5.10964563819
    iteration:  100 , logli:  5.11265838512
    Num of iterations:  107 , logli:  5.1126810697
    EM algorithm terminated. [1, 1, 1, 2]
    Trying orders  [1, 1, 1, 1, 1]
    iteration:  10 , logli:  5.03154595258
    iteration:  20 , logli:  5.04393763194
    Num of iterations:  23 , logli:  5.04394031758
    EM algorithm terminated. [1, 1, 1, 1, 1]
    Best solution: logli = 5.12395778503 orders = [1, 1, 3]
    >>> print(D0)
    [[  -83.42943     0.          0.          0.          0.     ]
     [    0.       -718.67799     0.          0.          0.     ]
     [    0.          0.      -1026.16063  1026.16063     0.     ]
     [    0.          0.          0.      -1026.16063  1026.16063]
     [    0.          0.          0.          0.      -1026.16063]]
    >>> print(D1)
    [[  54.14857    4.90186   24.379      0.         0.     ]
     [   3.39152  665.84735   49.43912    0.         0.     ]
     [   0.         0.         0.         0.         0.     ]
     [   0.         0.         0.         0.         0.     ]
     [  42.6473    96.94396  886.56937    0.         0.     ]]
    >>> logli = LikelihoodFromTrace(tr, D0, D1)
    >>> print(logli)
    5.123958173279035
    >>> trAcf = LagCorrelationsFromTrace(tr, 10)
    >>> print(trAcf)
    [0.18413533155701942, 0.18160338347883728, 0.17545482361568204, 0.19966682687727969, 0.083236073275010994, 0.086348372442235991, 0.095009981958434644, 0.062859937512461148, 0.062326277833923117, 0.065928707515509583]
    >>> mapAcf = LagCorrelationsFromMAP(D0, D1, 10)
    >>> print(mapAcf)
    [ 0.24889  0.17665  0.12882  0.09638  0.0738   0.05765  0.04578  0.03684  0.02995  0.02454]
    >>> sqAcf = SquaredDifference(mapAcf, trAcf)
    >>> print(sqAcf)
    0.0238340444119
    >>> reAcf = RelativeEntropy(mapAcf, trAcf)
    >>> print(reAcf)
    0.321362238531

