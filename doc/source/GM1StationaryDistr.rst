butools.mam.GM1StationaryDistr
==============================

.. currentmodule:: butools.mam

.. np:function:: GM1StationaryDistr

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`pi = GM1StationaryDistr (B, R, K)`
        * - Mathematica:
          - :code:`pi = GM1StationaryDistr [B, R, K]`
        * - Python/Numpy:
          - :code:`pi = GM1StationaryDistr (B, R, K)`

    Returns the stationary distribution of the G/M/1 type
    Markov chain up to a given level K.
    
    Parameters
    ----------
    A : length(M) list of matrices of shape (N,N)
        Matrix blocks of the G/M/1 type generator in the 
        regular part, from 0 to M-1.
    B : length(M) list of matrices of shape (N,N)
        Matrix blocks of the G/M/1 type generator at the
    R : matrix, shape (N,N)
        Matrix R of the G/M/1 type Markov chain
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

    >>> B0 = [0.7, 0.2; 0.3, 0.6];
    >>> B1 = [0.3, 0.4; 0.5, 0.2];
    >>> B2 = [0.2, 0.4; 0.1, 0.6];
    >>> B3 = [0., 0.1; 0.2, 0.];
    >>> A0 = [0.1, 0.; 0., 0.1];
    >>> A1 = [0., 0.2; 0., 0.2];
    >>> A2 = [0., 0.1; 0., 0.];
    >>> A3 = [0.3, 0.2; 0.3, 0.2];
    >>> A4 = [0., 0.1; 0.2, 0.];
    >>> B = {B0, B1, B2, B3};
    >>> A = {A0, A1, A2, A3, A4};
    >>> R = GM1FundamentalMatrix(A);
    >>> disp(R);
          0.10065     0.026961
       0.00065531      0.12569
    >>> pi = GM1StationaryDistr(B, R, 300);

    For Mathematica:

    >>> B0 = {{0.7, 0.2},{0.3, 0.6}};
    >>> B1 = {{0.3, 0.4},{0.5, 0.2}};
    >>> B2 = {{0.2, 0.4},{0.1, 0.6}};
    >>> B3 = {{0., 0.1},{0.2, 0.}};
    >>> A0 = {{0.1, 0.},{0., 0.1}};
    >>> A1 = {{0., 0.2},{0., 0.2}};
    >>> A2 = {{0., 0.1},{0., 0.}};
    >>> A3 = {{0.3, 0.2},{0.3, 0.2}};
    >>> A4 = {{0., 0.1},{0.2, 0.}};
    >>> B = {B0, B1, B2, B3};
    >>> A = {A0, A1, A2, A3, A4};
    >>> R = GM1FundamentalMatrix[A];
    "The evaluation of the iteration required "64" roots\n"
    "The evaluation of the iteration required "32" roots\n"
    "The evaluation of the iteration required "16" roots\n"
    "The evaluation of the iteration required "8" roots\n"
    "The evaluation of the iteration required "8" roots\n"
    "Final Residual Error for G: "5.551115123125783*^-17
    >>> Print[R];
    {{0.10065149910973312, 0.026960920607274754},
     {0.0006553100576153258, 0.12568710472819553}}
    >>> pi = GM1StationaryDistr[B, R, 300];
    "Accumulated mass after "2" iterations: "0.9838720044873233
    "Accumulated mass after "3" iterations: "0.9979548824322513
    "Accumulated mass after "4" iterations: "0.9997408547470504
    "Accumulated mass after "5" iterations: "0.9999671812477241
    "Accumulated mass after "6" iterations: "0.9999958456126867
    "Accumulated mass after "7" iterations: "0.999999474298702
    "Accumulated mass after "8" iterations: "0.9999999334955769
    "Accumulated mass after "9" iterations: "0.9999999915886283
    "Accumulated mass after "10" iterations: "0.9999999989363275
    "Accumulated mass after "11" iterations: "0.9999999998655101
    "Accumulated mass after "12" iterations: "0.999999999982997

    For Python/Numpy:

    >>> B0 = ml.matrix([[0.7, 0.2],[0.3, 0.6]])
    >>> B1 = ml.matrix([[0.3, 0.4],[0.5, 0.2]])
    >>> B2 = ml.matrix([[0.2, 0.4],[0.1, 0.6]])
    >>> B3 = ml.matrix([[0., 0.1],[0.2, 0.]])
    >>> A0 = ml.matrix([[0.1, 0.],[0., 0.1]])
    >>> A1 = ml.matrix([[0., 0.2],[0., 0.2]])
    >>> A2 = ml.matrix([[0., 0.1],[0., 0.]])
    >>> A3 = ml.matrix([[0.3, 0.2],[0.3, 0.2]])
    >>> A4 = ml.matrix([[0., 0.1],[0.2, 0.]])
    >>> B = [B0, B1, B2, B3]
    >>> A = [A0, A1, A2, A3, A4]
    >>> R = GM1FundamentalMatrix(A)
    The Shifted PWCR evaluation of Iteration  1  required  64  roots
    The Shifted PWCR evaluation of Iteration  2  required  32  roots
    The Shifted PWCR evaluation of Iteration  3  required  16  roots
    The Shifted PWCR evaluation of Iteration  4  required  8  roots
    The Shifted PWCR evaluation of Iteration  5  required  8  roots
    Final Residual Error for G:  5.20417042793e-17
    >>> print(R)
    [[ 0.10065  0.02696]
     [ 0.00066  0.12569]]
    >>> pi = GM1StationaryDistr(B, R, 300)
    Accumulated mass after  2  iterations:  0.983872004487
    Accumulated mass after  3  iterations:  0.997954882432
    Accumulated mass after  4  iterations:  0.999740854747
    Accumulated mass after  5  iterations:  0.999967181248
    Accumulated mass after  6  iterations:  0.999995845613
    Accumulated mass after  7  iterations:  0.999999474299
    Accumulated mass after  8  iterations:  0.999999933496
    Accumulated mass after  9  iterations:  0.999999991589
    Accumulated mass after  10  iterations:  0.999999998936
    Accumulated mass after  11  iterations:  0.999999999866
    Accumulated mass after  12  iterations:  0.999999999983

