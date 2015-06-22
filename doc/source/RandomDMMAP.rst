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

    >>> D = RandomDMMAP(4,3,5.62,10);
    >>> disp(D{1});
          0.56425            0    0.0025792     0.025159
        0.0054367      0.85152    0.0013472     0.010131
         0.013513    0.0018426      0.88354     0.013932
                0            0     0.032424      0.61785
    >>> disp(D{2});
         0.028968     0.028484     0.042247     0.033735
         0.013261     0.011258            0     0.011753
        0.0054317   0.00037558    0.0073518    0.0066517
         0.055321            0     0.027431     0.034022
    >>> disp(D{3});
         0.011736     0.039391     0.041982     0.030391
          0.01497    0.0067164     0.017958     0.012237
                0    0.0013491     0.015416    0.0043155
         0.027236            0     0.063043     0.033505
    >>> disp(D{4});
                0     0.061636     0.028295     0.061148
         0.018395            0    0.0081089     0.016913
         0.013366            0     0.015417     0.017502
         0.022411     0.014965     0.022161     0.049628
    >>> m = MarginalMomentsFromDMMAP(D,1);
    >>> disp(m);
             5.62

    For Mathematica:

    >>> D = RandomDMMAP[4,3,5.62,10];
    >>> Print[D[[1]]];
    {{0.7810890088422581, 0.018939409239872213, 0.011566728894108407, 0.012840905332734302},
     {0., 0.7686108676453307, 0., 0.0296274461881119},
     {0.006909444338897623, 0., 0.8212382065840111, 0.009918502382455511},
     {0.006725380233543529, 0.013446628384419759, 0.03029896637954113, 0.7411223072381175}}
    >>> Print[D[[2]]];
    {{0.021786721693507916, 0.001714487642362968, 0.03205664006843156, 0.03960954820224747},
     {0.02723766414211223, 0.01282358443867364, 0.0415239546381109, 0.03493927546794094},
     {0.009997402396806642, 0.005462615290309531, 0.024432102117094734, 0.010530689426270206},
     {0.02441240110092562, 0.00949641405494263, 0.016365908921523384, 0.0011060112138241045}}
    >>> Print[D[[3]]];
    {{0.0069238682446880005, 0., 0.011712615171458458, 0.},
     {0., 0., 0.03806938082210104, 0.0024381956968668054},
     {0.0229141930014625, 0.015036606419348106, 0.017357867464892204, 0.002624812177169526},
     {0.022630713471452606, 0.030801496231725692, 0.030330584575548266, 0.0058371788889074044}}
    >>> Print[D[[4]]];
    {{0., 0.013891383392792696, 0.044261745911886356, 0.0036069373636515626},
     {0.016091088219161945, 0.019071554302603246, 0.009566988438986684, 0.},
     {0.015849307687609455, 0.006343341806002841, 0.028925213600800453, 0.0024596953068695697},
     {0.02107558892301927, 0.02501344412263133, 0.02133697625987782, 0.}}
    >>> m = MarginalMomentsFromDMMAP[D,1][[1]];
    >>> Print[m];
    5.620000000000001

    For Python/Numpy:

    >>> D = RandomDMMAP(4,3,5.62,10)
    >>> print(D[0])
    [[  8.62309e-01   1.39367e-02   1.69215e-02   1.36871e-02]
     [  8.12661e-03   8.16361e-01   1.70023e-02   1.86885e-02]
     [  6.87350e-02   2.32146e-02   4.52912e-01   0.00000e+00]
     [  7.48325e-05   1.69236e-02   0.00000e+00   6.79616e-01]]
    >>> print(D[1])
    [[ 0.01449  0.00897  0.       0.     ]
     [ 0.0008   0.01882  0.       0.00527]
     [ 0.02305  0.06771  0.05505  0.06285]
     [ 0.01182  0.       0.       0.02009]]
    >>> print(D[2])
    [[ 0.00859  0.01967  0.00193  0.01385]
     [ 0.01537  0.00505  0.01802  0.01952]
     [ 0.04101  0.03806  0.03466  0.     ]
     [ 0.01504  0.018    0.0566   0.02218]]
    >>> print(D[3])
    [[ 0.00393  0.00258  0.01685  0.0023 ]
     [ 0.01851  0.02109  0.       0.01737]
     [ 0.05348  0.       0.06406  0.0152 ]
     [ 0.0667   0.03165  0.00101  0.0603 ]]
    >>> m = MarginalMomentsFromDMMAP(D,1)[0]
    >>> print(m)
    5.62

