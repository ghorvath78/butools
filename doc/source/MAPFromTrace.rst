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

