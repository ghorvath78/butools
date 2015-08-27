butools.fitting.PHFromTrace
===========================

.. currentmodule:: butools.fitting

.. np:function:: PHFromTrace

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A, logli] = PHFromTrace(trace, orders, maxIter, stopCond, initial, result)`
        * - Mathematica:
          - :code:`{alpha, A, logli} = PHFromTrace[trace, orders, maxIter, stopCond, initial, result]`
        * - Python/Numpy:
          - :code:`alpha, A, logli = PHFromTrace(trace, orders, maxIter, stopCond, initial, result)`

    Performs PH distribution fitting using the EM algorithm
    (G-FIT, [1]_).
    
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
    initial : tuple of two vectors, optional
        The initial values of the branch probabilities and
        rate parameters is given by this tuple. If not 
        given, a default initial guess is determined and 
        the algorithm starts from there.
    result : {"vecmat", "vecvec"}, optional
        The result can be returned two ways. If "vecmat" is
        selected, the result is returned in the classical
        representation of phase-type distributions, thus the
        initial vector and the generator matrix. 
        If "vecvec" is selected, two vectors are returned, 
        one holds the branch probabilities, and the second
        holds the rate parameters of the Erlang branches.
        The default value is "vecmat"

    Returns
    -------
    (alpha, A) : tuple of matrix, shape (1,M) and matrix, shape (M,M)
        If the "vecmat" result format is chosen, the function
        returns the initial probability vector and the
        generator matrix of the phase type distribution.
    (pi, lambda) : tuple of vector, length N and vector, length N
        If the "vecvec" result format is chosen, the function
        returns the vector of branch probabilities and the
        vector of branch rates in a tuple.
    logli : double
        The log-likelihood divided by the trace length
        
    Notes
    -----
    This procedure is quite fast in the supported 
    mathematical frameworks. If the maximum speed is
    needed, please use the multi-core optimized c++
    implementation called SPEM-FIT_.

    .. _SPEM-FIT: https://bitbucket.org/ghorvath78/spemfit

    References
    ----------
    .. [1] Thummler, Axel, Peter Buchholz, and Miklós Telek.
           A novel approach for fitting probability 
           distributions to real trace data with the EM 
           algorithm. Dependable Systems and Networks, 2005.

    Examples
    ========
    For Matlab:

    >>> tr = dlmread('/home/gabor/github/butools/test/data/bctrace.iat');
    >>> [alpha, A] = PHFromTrace(tr, 5);
    Num of iterations: 26, logli: 4.80924
    Num of iterations: 26, logli: 4.80924
    EM algorithm terminated. (orders=1,4)
    Num of iterations: 38, logli: 4.70421
    EM algorithm terminated. (orders=2,3)
    Num of iterations: 30, logli: 4.89844
    Num of iterations: 60, logli: 4.8989
    Num of iterations: 91, logli: 4.89925
    Num of iterations: 121, logli: 4.8996
    Num of iterations: 150, logli: 4.9
    Num of iterations: 179, logli: 4.90046
    Num of iterations: 201, logli: 4.90074
    EM algorithm terminated. (orders=1,1,3)
    Num of iterations: 30, logli: 4.85155
    Num of iterations: 58, logli: 4.85184
    Num of iterations: 87, logli: 4.85214
    Num of iterations: 119, logli: 4.85259
    Num of iterations: 150, logli: 4.85461
    Num of iterations: 179, logli: 4.88627
    Num of iterations: 192, logli: 4.91508
    EM algorithm terminated. (orders=1,2,2)
    Num of iterations: 22, logli: 4.88641
    Num of iterations: 43, logli: 4.88777
    Num of iterations: 64, logli: 4.88814
    Num of iterations: 84, logli: 4.88827
    Num of iterations: 104, logli: 4.88833
    Num of iterations: 124, logli: 4.88836
    Num of iterations: 145, logli: 4.88837
    Num of iterations: 153, logli: 4.88838
    EM algorithm terminated. (orders=1,1,1,2)
    Num of iterations: 17, logli: 4.85091
    Num of iterations: 27, logli: 4.85108
    EM algorithm terminated. (orders=1,1,1,1,1)
    Best solution: logli=4.91508, orders=1,2,2
    >>> disp(alpha);
         0.065027      0.85788            0     0.077088            0
    >>> disp(A);
          -63.308            0            0            0            0
                0      -815.72       815.72            0            0
                0            0      -815.72            0            0
                0            0            0       -12563        12563
                0            0            0            0       -12563
    >>> logli = LikelihoodFromTrace(tr, alpha, A);
    >>> disp(logli);
           4.9151
    >>> intBounds = linspace(0, MarginalMomentsFromTrace(tr, 1)*4, 50);
    >>> [pdfTrX, pdfTrY] = PdfFromTrace(tr, intBounds);
    >>> [pdfPHX, pdfPHY] = IntervalPdfFromPH(alpha, A, intBounds);
    >>> sqPdf = EmpiricalSquaredDifference(pdfTrY, pdfPHY, intBounds);
    >>> disp(sqPdf);
        0.0079115
    >>> rePdf = EmpiricalRelativeEntropy(pdfTrY, pdfPHY, intBounds);
    >>> disp(rePdf);
          0.35834
    >>> [cdfTrX, cdfTrY] = CdfFromTrace(tr);
    >>> step = ceil(length(tr)/2000);
    >>> cdfTrX = cdfTrX(1:step:length(tr));
    >>> cdfTrY = cdfTrY(1:step:length(tr));
    >>> cdfPHY = CdfFromPH(alpha, A, cdfTrX);
    >>> sqCdf = EmpiricalSquaredDifference(cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX);
    >>> disp(sqCdf);
       9.9902e-11
    >>> reCdf = EmpiricalRelativeEntropy(cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX);
    >>> disp(reCdf);
       0.00018583

    For Mathematica:

    >>> tr = Flatten[Import["/home/gabor/github/butools/test/data/bctrace.iat","CSV"]];
    >>> {alpha, A, logli} = PHFromTrace[tr, 5];
    "Num of iterations: "20", logli: "4.8092300067059535
    "Num of iterations: "26", logli: "4.809235435027485
    "EM algorithm terminated. (orders="{1, 4}")"
    "Num of iterations: "20", logli: "4.683317781477141
    "Num of iterations: "38", logli: "4.704214271672482
    "Num of iterations: "38", logli: "4.704214271672482
    "EM algorithm terminated. (orders="{2, 3}")"
    "Num of iterations: "13", logli: "4.894318246490096
    "Num of iterations: "26", logli: "4.898349065288064
    "Num of iterations: "39", logli: "4.8985999642994065
    "Num of iterations: "51", logli: "4.89878133343547
    "Num of iterations: "64", logli: "4.8989479972108
    "Num of iterations: "77", logli: "4.8990977181707995
    "Num of iterations: "90", logli: "4.89924035598373
    "Num of iterations: "102", logli: "4.899372417693129
    "Num of iterations: "114", logli: "4.899510325011154
    "Num of iterations: "127", logli: "4.899671843230089
    "Num of iterations: "139", logli: "4.899836025597829
    "Num of iterations: "152", logli: "4.900031245080798
    "Num of iterations: "165", logli: "4.900238803829762
    "Num of iterations: "178", logli: "4.900443966924616
    "Num of iterations: "190", logli: "4.9006139623651945
    "Num of iterations: "201", logli: "4.900742157090518
    "EM algorithm terminated. (orders="{1, 1, 3}")"
    "Num of iterations: "14", logli: "4.846111351537217
    "Num of iterations: "27", logli: "4.851422822109602
    "Num of iterations: "40", logli: "4.851663586226162
    "Num of iterations: "52", logli: "4.85177202585639
    "Num of iterations: "65", logli: "4.851912954098613
    "Num of iterations: "78", logli: "4.8520519001971145
    "Num of iterations: "91", logli: "4.852183146393526
    "Num of iterations: "104", logli: "4.852333485343567
    "Num of iterations: "117", logli: "4.852550797457522
    "Num of iterations: "130", logli: "4.852926291581028
    "Num of iterations: "143", logli: "4.853720321486975
    "Num of iterations: "158", logli: "4.85667738204316
    "Num of iterations: "167", logli: "4.8623051325434705
    "Num of iterations: "168", logli: "4.863279856542737
    "Num of iterations: "169", logli: "4.864341586890624
    "Num of iterations: "170", logli: "4.865497327898647
    "Num of iterations: "171", logli: "4.866757372059696
    "Num of iterations: "172", logli: "4.868137493672351
    "Num of iterations: "173", logli: "4.8696621585039255
    "Num of iterations: "174", logli: "4.871369345657646
    "Num of iterations: "175", logli: "4.8733180347788565
    "Num of iterations: "176", logli: "4.875600054271152
    "Num of iterations: "177", logli: "4.878358389000584
    "Num of iterations: "178", logli: "4.881811700896143
    "Num of iterations: "179", logli: "4.886268931740982
    "Num of iterations: "180", logli: "4.892054411447398
    "Num of iterations: "181", logli: "4.899120695311693
    "Num of iterations: "182", logli: "4.9062691521540565
    "Num of iterations: "183", logli: "4.911394072570149
    "Num of iterations: "184", logli: "4.913841789490561
    "Num of iterations: "185", logli: "4.914707340661588
    "Num of iterations: "186", logli: "4.914972473547811
    "Num of iterations: "187", logli: "4.915049918282578
    "Num of iterations: "188", logli: "4.915072471899828
    "Num of iterations: "189", logli: "4.915079202008659
    "Num of iterations: "190", logli: "4.915081312881225
    "Num of iterations: "191", logli: "4.915082026039371
    "Num of iterations: "192", logli: "4.91508228990118
    "Num of iterations: "192", logli: "4.91508228990118
    "EM algorithm terminated. (orders="{1, 2, 2}")"
    "Num of iterations: "11", logli: "4.879855265152097
    "Num of iterations: "21", logli: "4.886270915069811
    "Num of iterations: "30", logli: "4.887168456960867
    "Num of iterations: "40", logli: "4.887668315460644
    "Num of iterations: "49", logli: "4.887917322768672
    "Num of iterations: "59", logli: "4.888081691124073
    "Num of iterations: "69", logli: "4.888181621001137
    "Num of iterations: "78", logli: "4.888239956657024
    "Num of iterations: "88", logli: "4.888283611362753
    "Num of iterations: "98", logli: "4.888313275634503
    "Num of iterations: "108", logli: "4.888334008025995
    "Num of iterations: "118", logli: "4.888348859128792
    "Num of iterations: "128", logli: "4.888359731068756
    "Num of iterations: "138", logli: "4.888367845219129
    "Num of iterations: "149", logli: "4.888374537830832
    "Num of iterations: "153", logli: "4.8883765319279915
    "EM algorithm terminated. (orders="{1, 1, 1, 2}")"
    "Num of iterations: "9", logli: "4.844577387233243
    "Num of iterations: "17", logli: "4.850912188288292
    "Num of iterations: "26", logli: "4.8510832077549315
    "Num of iterations: "27", logli: "4.851083604173109
    "EM algorithm terminated. (orders="{1, 1, 1, 1, 1}")"
    "Best solution: logli="4.91508228990118", orders="{1, 2, 2}
    >>> Print[alpha];
    {0.06502731323053666, 0.8578848706972437, 0, 0.07708781607263537, 0}
    >>> Print[A];
    {{-63.307631149011584, 0, 0, 0, 0},
     {0, -815.7180689395599, 815.7180689395599, 0, 0},
     {0, 0., -815.7180689395599, 0, 0},
     {0, 0, 0, -12563.087922793351, 12563.087922793351},
     {0, 0, 0, 0., -12563.087922793351}}
    >>> Print[logli];
    4.91508228990118
    >>> logli = LikelihoodFromTrace[tr, alpha, A];
    >>> Print[logli];
    4.915082396859105
    >>> intBounds = Array[# &, 50, {0, MarginalMomentsFromTrace[tr, 1][[1]]*4}];
    >>> {pdfTrX, pdfTrY} = PdfFromTrace[tr, intBounds];
    >>> {pdfPHX, pdfPHY} = IntervalPdfFromPH[alpha, A, intBounds];
    >>> sqPdf = EmpiricalSquaredDifference[pdfTrY, pdfPHY, intBounds];
    >>> Print[sqPdf];
    0.007911510224685201
    >>> rePdf = EmpiricalRelativeEntropy[pdfTrY, pdfPHY, intBounds];
    >>> Print[rePdf];
    0.3583352563975375
    >>> {cdfTrX, cdfTrY} = CdfFromTrace[tr];
    >>> step = Ceiling[Length[tr]/2000];
    >>> cdfTrX = cdfTrX[[1;;Length[tr];;step]];
    >>> cdfTrY = cdfTrY[[1;;Length[tr];;step]];
    >>> cdfPHY = CdfFromPH[alpha, A, cdfTrX];
    >>> sqCdf = EmpiricalSquaredDifference[cdfTrY[[1;;-2]], cdfPHY[[1;;-2]], cdfTrX];
    >>> Print[sqCdf];
    9.990199283948094*^-11
    >>> reCdf = EmpiricalRelativeEntropy[cdfTrY[[1;;-2]], cdfPHY[[1;;-2]], cdfTrX];
    >>> Print[reCdf];
    0.000185826390423264

    For Python/Numpy:

    >>> tr = np.loadtxt("/home/gabor/github/butools/test/data/bctrace.iat")
    >>> alpha, A = PHFromTrace(tr, 5)
    iteration:  10 , logli:  4.8091271432
    iteration:  20 , logli:  4.80923172357
    EM algorithm terminated. [1, 4]
    Num of iterations:  26 , logli:  4.80923543503
    iteration:  10 , logli:  4.65521504335
    iteration:  20 , logli:  4.68768389352
    iteration:  30 , logli:  4.70415249872
    EM algorithm terminated. [2, 3]
    Num of iterations:  38 , logli:  4.70421427167
    iteration:  10 , logli:  4.88937597779
    iteration:  20 , logli:  4.89819514171
    iteration:  30 , logli:  4.89845514506
    iteration:  40 , logli:  4.89863272893
    iteration:  50 , logli:  4.89878133343
    iteration:  60 , logli:  4.89891146874
    iteration:  70 , logli:  4.89903002265
    iteration:  80 , logli:  4.89914197577
    iteration:  90 , logli:  4.89925126394
    iteration:  100 , logli:  4.89936125929
    iteration:  110 , logli:  4.89947502785
    iteration:  120 , logli:  4.89959540205
    iteration:  130 , logli:  4.89972483134
    iteration:  140 , logli:  4.89986492887
    iteration:  150 , logli:  4.90001566233
    iteration:  160 , logli:  4.90017434189
    iteration:  170 , logli:  4.90033497306
    iteration:  180 , logli:  4.90048885472
    iteration:  190 , logli:  4.9006268209
    iteration:  200 , logli:  4.90074215709
    EM algorithm terminated. [1, 1, 3]
    Num of iterations:  201 , logli:  4.90074215709
    iteration:  10 , logli:  4.8428939533
    iteration:  20 , logli:  4.85044993659
    iteration:  30 , logli:  4.85156877825
    iteration:  40 , logli:  4.85167161239
    iteration:  50 , logli:  4.85176192227
    iteration:  60 , logli:  4.8518686018
    iteration:  70 , logli:  4.85197855253
    iteration:  80 , logli:  4.85208231474
    iteration:  90 , logli:  4.8521831464
    iteration:  100 , logli:  4.85229494225
    iteration:  110 , logli:  4.85243806034
    iteration:  120 , logli:  4.85264283694
    iteration:  130 , logli:  4.85296689978
    iteration:  140 , logli:  4.85354745937
    iteration:  150 , logli:  4.85478414915
    iteration:  160 , logli:  4.85802067972
    iteration:  170 , logli:  4.86675737206
    iteration:  180 , logli:  4.89912069533
    iteration:  190 , logli:  4.91508202604
    EM algorithm terminated. [1, 2, 2]
    Num of iterations:  192 , logli:  4.9150822899
    iteration:  10 , logli:  4.87985526515
    iteration:  20 , logli:  4.88627091507
    iteration:  30 , logli:  4.88723479659
    iteration:  40 , logli:  4.88770291202
    iteration:  50 , logli:  4.88795746078
    iteration:  60 , logli:  4.88810565844
    iteration:  70 , logli:  4.88819663412
    iteration:  80 , logli:  4.8882549448
    iteration:  90 , logli:  4.88829369783
    iteration:  100 , logli:  4.88832026537
    iteration:  110 , logli:  4.88833897739
    iteration:  120 , logli:  4.88835247272
    iteration:  130 , logli:  4.8883624119
    iteration:  140 , logli:  4.88836986989
    iteration:  150 , logli:  4.88837556031
    EM algorithm terminated. [1, 1, 1, 2]
    Num of iterations:  153 , logli:  4.88837653192
    iteration:  10 , logli:  4.84789680934
    iteration:  20 , logli:  4.85106678625
    EM algorithm terminated. [1, 1, 1, 1, 1]
    Num of iterations:  27 , logli:  4.85108360418
    >>> print(alpha)
    [[ 0.06503  0.85788  0.       0.07709  0.     ]]
    >>> print(A)
    [[   -63.30763      0.           0.           0.           0.     ]
     [     0.        -815.71807    815.71807      0.           0.     ]
     [     0.           0.        -815.71807      0.           0.     ]
     [     0.           0.           0.      -12563.08792  12563.08792]
     [     0.           0.           0.           0.      -12563.08792]]
    >>> logli = LikelihoodFromTrace(tr, alpha, A)
    >>> print(logli)
    4.91508239686
    >>> intBounds = np.linspace(0, MarginalMomentsFromTrace(tr, 1)[0]*4, 50)
    >>> pdfTrX, pdfTrY = PdfFromTrace(tr, intBounds)
    >>> pdfPHX, pdfPHY = IntervalPdfFromPH(alpha, A, intBounds)
    >>> sqPdf = EmpiricalSquaredDifference(pdfTrY, pdfPHY, intBounds)
    >>> print(sqPdf)
    0.00791151022468
    >>> rePdf = EmpiricalRelativeEntropy(pdfTrY, pdfPHY, intBounds)
    >>> print(rePdf)
    0.358335256398
    >>> cdfTrX, cdfTrY = CdfFromTrace(tr)
    >>> step = math.ceil(Length(tr)/2000)
    >>> cdfTrX = cdfTrX[0:Length(tr):step]
    >>> cdfTrY = cdfTrY[0:Length(tr):step]
    >>> cdfPHY = CdfFromPH(alpha, A, cdfTrX)
    >>> sqCdf = EmpiricalSquaredDifference(cdfTrY[0:-1], cdfPHY[0:-1], cdfTrX)
    >>> print(sqCdf)
    9.99019928404e-11
    >>> reCdf = EmpiricalRelativeEntropy(cdfTrY[0:-1], cdfPHY[0:-1], cdfTrX)
    >>> print(reCdf)
    0.000185826390424

