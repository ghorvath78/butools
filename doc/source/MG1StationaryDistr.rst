butools.mam.MG1StationaryDistr
==============================

.. currentmodule:: butools.mam

.. np:function:: MG1StationaryDistr

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pi = MG1StationaryDistr (A, B, G, K)`
        * - Mathematica:
          - :code:`pi = MG1StationaryDistr [A, B, G, K]`
        * - Python/Numpy:
          - :code:`pi = MG1StationaryDistr (A, B, G, K)`

    Returns the stationary distribution of the M/G/1 type
    Markov chain up to a given level K.
    
    Parameters
    ----------
    A : length(M) list of matrices of shape (N,N)
        Matrix blocks of the M/G/1 type generator in the 
        regular part, from 0 to M-1.
    B : length(M) list of matrices of shape (N,N)
        Matrix blocks of the M/G/1 type generator at the
        boundary, from 0 to M-1.
    G : matrix, shape (N,N)
        Matrix G of the M/G/1 type Markov chain
    K : integer
        The stationary distribution is returned up to
        this level.
    
    Returns
    -------
    pi : array, shape (1,(K+1)*N)
        The stationary probability vector up to level K
    
    Examples
    ========
    For Matlab:

    >>> B0 = [0.1, 0.5; 0.3, 0.4];
    >>> B1 = [0., 0.1; 0., 0.];
    >>> B2 = [0.2, 0.; 0., 0.2];
    >>> B3 = [0., 0.1; 0.1, 0.];
    >>> A0 = [0.4, 0.2; 0.3, 0.4];
    >>> A1 = [0., 0.1; 0., 0.];
    >>> A2 = [0., 0.2; 0., 0.2];
    >>> A3 = [0.1, 0.; 0.1, 0.];
    >>> B = {B0, B1, B2, B3};
    >>> A = {A0, A1, A2, A3};
    >>> G = MG1FundamentalMatrix(A);
    >>> disp(G);
          0.60503      0.39497
          0.45912      0.54088
    >>> pi = MG1StationaryDistr(A, B, G, 300);

    For Mathematica:

    >>> B0 = {{0.1, 0.5},{0.3, 0.4}};
    >>> B1 = {{0., 0.1},{0., 0.}};
    >>> B2 = {{0.2, 0.},{0., 0.2}};
    >>> B3 = {{0., 0.1},{0.1, 0.}};
    >>> A0 = {{0.4, 0.2},{0.3, 0.4}};
    >>> A1 = {{0., 0.1},{0., 0.}};
    >>> A2 = {{0., 0.2},{0., 0.2}};
    >>> A3 = {{0.1, 0.},{0.1, 0.}};
    >>> B = {B0, B1, B2, B3};
    >>> A = {A0, A1, A2, A3};
    >>> G = MG1FundamentalMatrix[A];
    "The evaluation of the iteration required "64" roots\n"
    "The evaluation of the iteration required "32" roots\n"
    "The evaluation of the iteration required "16" roots\n"
    "The evaluation of the iteration required "16" roots\n"
    "The evaluation of the iteration required "8" roots\n"
    "Final Residual Error for G: "1.6653345369377348*^-16
    >>> Print[G];
    {{0.6050345283244288, 0.39496547167557117},
     {0.4591222984767535, 0.5408777015232465}}
    >>> pi = MG1StationaryDistr[A, B, G, 300];
    "Accumulated mass of the first "1" (reblocked) components:"0.3916489588121354
    "Accumulated mass of the first "2" (reblocked) components:"0.5693301777261455
    "Accumulated mass of the first "3" (reblocked) components:"0.7100355384481125
    "Accumulated mass of the first "4" (reblocked) components:"0.80138806057205
    "Accumulated mass of the first "5" (reblocked) components:"0.8647665851776781
    "Accumulated mass of the first "6" (reblocked) components:"0.9077316038451807
    "Accumulated mass of the first "7" (reblocked) components:"0.9370907444252965
    "Accumulated mass of the first "8" (reblocked) components:"0.9570975394581935
    "Accumulated mass of the first "9" (reblocked) components:"0.9707441149817535
    "Accumulated mass of the first "10" (reblocked) components:"0.9800493560702122
    "Accumulated mass of the first "11" (reblocked) components:"0.9863950724614493
    "Accumulated mass of the first "12" (reblocked) components:"0.9907223698850346
    "Accumulated mass of the first "13" (reblocked) components:"0.9936732983453151
    "Accumulated mass of the first "14" (reblocked) components:"0.9956856255466052
    "Accumulated mass of the first "15" (reblocked) components:"0.9970578944152654
    "Accumulated mass of the first "16" (reblocked) components:"0.9979936869664346
    "Accumulated mass of the first "17" (reblocked) components:"0.9986318329494369
    "Accumulated mass of the first "18" (reblocked) components:"0.9990670044714363
    "Accumulated mass of the first "19" (reblocked) components:"0.9993637614250387
    "Accumulated mass of the first "20" (reblocked) components:"0.9995661291912358
    "Accumulated mass of the first "21" (reblocked) components:"0.9997041300448283
    "Accumulated mass of the first "22" (reblocked) components:"0.9997982371051191
    "Accumulated mass of the first "23" (reblocked) components:"0.9998624116270072
    "Accumulated mass of the first "24" (reblocked) components:"0.9999061742229949
    "Accumulated mass of the first "25" (reblocked) components:"0.9999360172939101
    "Accumulated mass of the first "26" (reblocked) components:"0.9999563682091505
    "Accumulated mass of the first "27" (reblocked) components:"0.999970246129164
    "Accumulated mass of the first "28" (reblocked) components:"0.9999797099130591
    "Accumulated mass of the first "29" (reblocked) components:"0.9999861635606897
    "Accumulated mass of the first "30" (reblocked) components:"0.9999905645030824
    "Accumulated mass of the first "31" (reblocked) components:"0.9999935656421364
    "Accumulated mass of the first "32" (reblocked) components:"0.9999956122118974
    "Accumulated mass of the first "33" (reblocked) components:"0.9999970078312643
    "Accumulated mass of the first "34" (reblocked) components:"0.9999979595473771
    "Accumulated mass of the first "35" (reblocked) components:"0.9999986085520992
    "Accumulated mass of the first "36" (reblocked) components:"0.9999990511285394
    "Accumulated mass of the first "37" (reblocked) components:"0.9999993529351345
    "Accumulated mass of the first "38" (reblocked) components:"0.9999995587464082
    "Accumulated mass of the first "39" (reblocked) components:"0.9999996990954962
    "Accumulated mass of the first "40" (reblocked) components:"0.9999997948038905
    "Accumulated mass of the first "41" (reblocked) components:"0.9999998600704115
    "Accumulated mass of the first "42" (reblocked) components:"0.9999999045776756
    "Accumulated mass of the first "43" (reblocked) components:"0.9999999349285588
    "Accumulated mass of the first "44" (reblocked) components:"0.9999999556257669
    "Accumulated mass of the first "45" (reblocked) components:"0.9999999697398348
    "Accumulated mass of the first "46" (reblocked) components:"0.9999999793646553
    "Accumulated mass of the first "47" (reblocked) components:"0.9999999859281188
    "Accumulated mass of the first "48" (reblocked) components:"0.999999990403948
    "Accumulated mass of the first "49" (reblocked) components:"0.9999999934561548
    "Accumulated mass of the first "50" (reblocked) components:"0.9999999955375491
    "Accumulated mass of the first "51" (reblocked) components:"0.9999999969569165
    "Accumulated mass of the first "52" (reblocked) components:"0.9999999979248272
    "Accumulated mass of the first "53" (reblocked) components:"0.9999999985848754
    "Accumulated mass of the first "54" (reblocked) components:"0.9999999990349828
    "Accumulated mass of the first "55" (reblocked) components:"0.999999999341925
    "Accumulated mass of the first "56" (reblocked) components:"0.9999999995512383
    "Accumulated mass of the first "57" (reblocked) components:"0.9999999996939756
    "Accumulated mass of the first "58" (reblocked) components:"0.9999999997913126
    "Accumulated mass of the first "59" (reblocked) components:"0.9999999998576897
    "Accumulated mass of the first "60" (reblocked) components:"0.9999999999029543

    For Python/Numpy:

    >>> B0 = ml.matrix([[0.1, 0.5],[0.3, 0.4]])
    >>> B1 = ml.matrix([[0., 0.1],[0., 0.]])
    >>> B2 = ml.matrix([[0.2, 0.],[0., 0.2]])
    >>> B3 = ml.matrix([[0., 0.1],[0.1, 0.]])
    >>> A0 = ml.matrix([[0.4, 0.2],[0.3, 0.4]])
    >>> A1 = ml.matrix([[0., 0.1],[0., 0.]])
    >>> A2 = ml.matrix([[0., 0.2],[0., 0.2]])
    >>> A3 = ml.matrix([[0.1, 0.],[0.1, 0.]])
    >>> B = [B0, B1, B2, B3]
    >>> A = [A0, A1, A2, A3]
    >>> G = MG1FundamentalMatrix(A)
    The Shifted PWCR evaluation of Iteration  1  required  64  roots
    The Shifted PWCR evaluation of Iteration  2  required  32  roots
    The Shifted PWCR evaluation of Iteration  3  required  16  roots
    The Shifted PWCR evaluation of Iteration  4  required  16  roots
    The Shifted PWCR evaluation of Iteration  5  required  8  roots
    Final Residual Error for G:  1.66533453694e-16
    >>> print(G)
    [[ 0.60503  0.39497]
     [ 0.45912  0.54088]]
    >>> pi = MG1StationaryDistr(A, B, G, 300)
    Accumulated mass of the first  1  (reblocked) components: 0.391648958812
    Accumulated mass of the first  2  (reblocked) components: 0.569330177726
    Accumulated mass of the first  3  (reblocked) components: 0.710035538448
    Accumulated mass of the first  4  (reblocked) components: 0.801388060572
    Accumulated mass of the first  5  (reblocked) components: 0.864766585178
    Accumulated mass of the first  6  (reblocked) components: 0.907731603845
    Accumulated mass of the first  7  (reblocked) components: 0.937090744425
    Accumulated mass of the first  8  (reblocked) components: 0.957097539458
    Accumulated mass of the first  9  (reblocked) components: 0.970744114982
    Accumulated mass of the first  10  (reblocked) components: 0.98004935607
    Accumulated mass of the first  11  (reblocked) components: 0.986395072461
    Accumulated mass of the first  12  (reblocked) components: 0.990722369885
    Accumulated mass of the first  13  (reblocked) components: 0.993673298345
    Accumulated mass of the first  14  (reblocked) components: 0.995685625547
    Accumulated mass of the first  15  (reblocked) components: 0.997057894415
    Accumulated mass of the first  16  (reblocked) components: 0.997993686966
    Accumulated mass of the first  17  (reblocked) components: 0.998631832949
    Accumulated mass of the first  18  (reblocked) components: 0.999067004471
    Accumulated mass of the first  19  (reblocked) components: 0.999363761425
    Accumulated mass of the first  20  (reblocked) components: 0.999566129191
    Accumulated mass of the first  21  (reblocked) components: 0.999704130045
    Accumulated mass of the first  22  (reblocked) components: 0.999798237105
    Accumulated mass of the first  23  (reblocked) components: 0.999862411627
    Accumulated mass of the first  24  (reblocked) components: 0.999906174223
    Accumulated mass of the first  25  (reblocked) components: 0.999936017294
    Accumulated mass of the first  26  (reblocked) components: 0.999956368209
    Accumulated mass of the first  27  (reblocked) components: 0.999970246129
    Accumulated mass of the first  28  (reblocked) components: 0.999979709913
    Accumulated mass of the first  29  (reblocked) components: 0.999986163561
    Accumulated mass of the first  30  (reblocked) components: 0.999990564503
    Accumulated mass of the first  31  (reblocked) components: 0.999993565642
    Accumulated mass of the first  32  (reblocked) components: 0.999995612212
    Accumulated mass of the first  33  (reblocked) components: 0.999997007831
    Accumulated mass of the first  34  (reblocked) components: 0.999997959547
    Accumulated mass of the first  35  (reblocked) components: 0.999998608552
    Accumulated mass of the first  36  (reblocked) components: 0.999999051129
    Accumulated mass of the first  37  (reblocked) components: 0.999999352935
    Accumulated mass of the first  38  (reblocked) components: 0.999999558746
    Accumulated mass of the first  39  (reblocked) components: 0.999999699095
    Accumulated mass of the first  40  (reblocked) components: 0.999999794804
    Accumulated mass of the first  41  (reblocked) components: 0.99999986007
    Accumulated mass of the first  42  (reblocked) components: 0.999999904578
    Accumulated mass of the first  43  (reblocked) components: 0.999999934929
    Accumulated mass of the first  44  (reblocked) components: 0.999999955626
    Accumulated mass of the first  45  (reblocked) components: 0.99999996974
    Accumulated mass of the first  46  (reblocked) components: 0.999999979365
    Accumulated mass of the first  47  (reblocked) components: 0.999999985928
    Accumulated mass of the first  48  (reblocked) components: 0.999999990404
    Accumulated mass of the first  49  (reblocked) components: 0.999999993456
    Accumulated mass of the first  50  (reblocked) components: 0.999999995538
    Accumulated mass of the first  51  (reblocked) components: 0.999999996957
    Accumulated mass of the first  52  (reblocked) components: 0.999999997925
    Accumulated mass of the first  53  (reblocked) components: 0.999999998585
    Accumulated mass of the first  54  (reblocked) components: 0.999999999035
    Accumulated mass of the first  55  (reblocked) components: 0.999999999342
    Accumulated mass of the first  56  (reblocked) components: 0.999999999551
    Accumulated mass of the first  57  (reblocked) components: 0.999999999694
    Accumulated mass of the first  58  (reblocked) components: 0.999999999791
    Accumulated mass of the first  59  (reblocked) components: 0.999999999858
    Accumulated mass of the first  60  (reblocked) components: 0.999999999903

