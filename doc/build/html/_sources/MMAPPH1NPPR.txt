butools.queues.MMAPPH1NPPR
==========================

.. currentmodule:: butools.queues

.. np:function:: MMAPPH1NPPR

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Ret = MMAPPH1NPPR(D, sigma, S, ...)`
        * - Mathematica:
          - :code:`Ret = MMAPPH1NPPR[D, sigma, S, ...]`
        * - Python/Numpy:
          - :code:`Ret = MMAPPH1NPPR(D, sigma, S, ...)`

    Returns various performane measures of a continuous time 
    MMAP[K]/PH[K]/1 non-preemptive priority queue, see [1]_.

    Parameters
    ----------
    D : list of matrices of shape (N,N), length (K+1)
        The D0...DK matrices of the arrival process.
        D1 corresponds to the lowest, DK to the highest priority.
    sigma : list of row vectors, length (K)
        The list containing the initial probability vectors of the service
        time distributions of the various customer types. The length of the
       vectors does not have to be the same.
    S : list of square matrices, length (K)
        The transient generators of the phase type distributions representing
        the service time of the jobs belonging to various types.
    further parameters : 
        The rest of the function parameters specify the options
        and the performance measures to be computed.

        The supported performance measures and options in this 
        function are:

        +----------------+--------------------+----------------------------------------+
        | Parameter name | Input parameters   | Output                                 |
        +================+====================+========================================+
        | "ncMoms"       | Number of moments  | The moments of the number of customers |
        +----------------+--------------------+----------------------------------------+
        | "ncDistr"      | Upper limit K      | The distribution of the number of      |
        |                |                    | customers from level 0 to level K-1    |
        +----------------+--------------------+----------------------------------------+
        | "stMoms"       | Number of moments  | The sojourn time moments               |
        +----------------+--------------------+----------------------------------------+
        | "stDistr"      | A vector of points | The sojourn time distribution at the   |
        |                |                    | requested points (cummulative, cdf)    |
        +----------------+--------------------+----------------------------------------+
        | "prec"         | The precision      | Numerical precision used as a stopping |
        |                |                    | condition when solving the Riccati and |
        |                |                    | the matrix-quadratic equations         |
        +----------------+--------------------+----------------------------------------+
        | "erlMaxOrder"  | Integer number     | The maximal Erlang order used in the   |
        |                |                    | erlangization procedure. The default   |
        |                |                    | value is 200.                          |
        +----------------+--------------------+----------------------------------------+
        | "classes"      | Vector of integers | Only the performance measures          |
        |                |                    | belonging to these classes are         |
        |                |                    | returned. If not given, all classes    |
        |                |                    | are analyzed.                          |
        +----------------+--------------------+----------------------------------------+
        
        (The quantities related to the number of customers in 
        the system include the customer in the server, and the 
        sojourn time related quantities include the service 
        times as well)

    Returns
    -------
    Ret : list of the performance measures
        Each entry of the list corresponds to a performance 
        measure requested. Each entry is a matrix, where the
        columns belong to the various job types.
        If there is just a single item, 
        then it is not put into a list.

    References
    ----------
    .. [1] G. Horvath, "Efficient analysis of the MMAP[K]/PH[K]/1
           priority queue", European Journal of Operational 
           Research, 246(1), 128-139, 2015.

    Examples
    ========
    For Matlab:

    >>> D0 = [-5.49, 0., 1.15, 0.; 0., -2.29, 0., 0.; 0., 0.08, -1.32, 0.; 0.72, 1.17, 0.7, -7.07];
    >>> D1 = [0.25, 0.38, 0.64, 0.; 0., 0., 0., 1.09; 0., 1.24, 0., 0.; 0.37, 0., 0., 0.];
    >>> D2 = [0.3, 1.0, 0., 0.48; 0., 0.2, 0., 0.; 0., 0., 0., 0.; 0.61, 0., 0., 0.2];
    >>> D3 = [0., 0.98, 0., 0.31; 0., 0., 1.0, 0.; 0., 0., 0., 0.; 1.1, 0.84, 0.33, 1.03];
    >>> sigma3 = [0.83333,0.11404,0.05263];
    >>> S3 = [-3., 0., 0.; 0.73077, -0.73077, 0.; 0., 0.5, -0.5];
    >>> sigma2 = [1.];
    >>> S2 = [-2.];
    >>> sigma1 = [0.25,0.75];
    >>> S1 = [-2.5, 2.5; 0., -10.];
    >>> [ncm1, ncd1, ncm2, ncd2, ncm3, ncd3] = MMAPPH1NPPR({D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, 'ncMoms', 3, 'ncDistr', 500);
    Final Residual Error for Psi:     4.774e-15
    Final Residual Error for Psi:     4.774e-15
    Final Residual Error for Psi:    2.7894e-15
    Final Residual Error for Psi:    3.5458e-15
    Final Residual Error for Psi:    6.5281e-14
    Final Residual Error for Psi:    2.5535e-15
    Final Residual Error for Psi:    7.2442e-15
    Final Residual Error for Psi:    4.8017e-15
    >>> distrPoints = [1., 5., 10.];
    >>> [stm1, std1, stm2, std2, stm3, std3] = MMAPPH1NPPR({D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, 'stMoms', 3, 'stDistr', distrPoints);
    Final Residual Error for Psi:     4.774e-15
    Final Residual Error for Psi:     4.774e-15
    Final Residual Error for Psi:    2.7894e-15
    Final Residual Error for Psi:    3.5458e-15
    Final Residual Error for Psi:    6.5281e-14
    Final Residual Error for Psi:    9.1731e-16
    Final Residual Error for Psi:    1.1933e-15
    Final Residual Error for Psi:    1.3175e-15
    Final Residual Error for Psi:    7.2442e-15
    Final Residual Error for Psi:    3.7603e-15
    Final Residual Error for Psi:    2.3189e-15
    Final Residual Error for Psi:    2.0864e-15
    >>> disp(stm1);
           15.909       788.48        63966
    >>> disp(std1);
          0.24787      0.44649      0.57919
    >>> disp(stm2);
           5.3735       102.75         3651
    >>> disp(std2);
          0.32134      0.70548      0.83911
    >>> disp(stm3);
           2.2552       13.114       124.73
    >>> disp(std3);
          0.45672      0.86989      0.97222

    For Mathematica:

    >>> D0 = {{-5.49, 0., 1.15, 0.},{0., -2.29, 0., 0.},{0., 0.08, -1.32, 0.},{0.72, 1.17, 0.7, -7.07}};
    >>> D1 = {{0.25, 0.38, 0.64, 0.},{0., 0., 0., 1.09},{0., 1.24, 0., 0.},{0.37, 0., 0., 0.}};
    >>> D2 = {{0.3, 1.0, 0., 0.48},{0., 0.2, 0., 0.},{0., 0., 0., 0.},{0.61, 0., 0., 0.2}};
    >>> D3 = {{0., 0.98, 0., 0.31},{0., 0., 1.0, 0.},{0., 0., 0., 0.},{1.1, 0.84, 0.33, 1.03}};
    >>> sigma3 = {0.83333,0.11404,0.05263};
    >>> S3 = {{-3., 0., 0.},{0.73077, -0.73077, 0.},{0., 0.5, -0.5}};
    >>> sigma2 = {1.};
    >>> S2 = {{-2.}};
    >>> sigma1 = {0.25,0.75};
    >>> S1 = {{-2.5, 2.5},{0., -10.}};
    >>> {ncm1, ncd1, ncm2, ncd2, ncm3, ncd3} = MMAPPH1NPPR[{D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, "ncMoms", 3, "ncDistr", 500];
    "Final Residual Error for Psi: "6.217248937900877*^-15
    "Final Residual Error for Psi: "6.217248937900877*^-15
    "Final Residual Error for Psi: "3.209238430557093*^-15
    "Final Residual Error for Psi: "3.587408148320037*^-15
    "Final Residual Error for Psi: "7.260858581048524*^-14
    "Final Residual Error for Psi: "3.0236230186275748*^-15
    "Final Residual Error for Psi: "4.7878367936959876*^-15
    "Final Residual Error for Psi: "4.884981308350689*^-15
    "Final Residual Error for G: "3.752092530320891*^-16
    "Final Residual Error for R: "4.110680098199597*^-16
    "Final Residual Error for G: "3.752092530320891*^-16
    "Final Residual Error for R: "4.110680098199597*^-16
    >>> distrPoints = {1., 5., 10.};
    >>> {stm1, std1, stm2, std2, stm3, std3} = MMAPPH1NPPR[{D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, "stMoms", 3, "stDistr", distrPoints];
    "Final Residual Error for Psi: "6.217248937900877*^-15
    "Final Residual Error for Psi: "6.217248937900877*^-15
    "Final Residual Error for Psi: "3.209238430557093*^-15
    "Final Residual Error for Psi: "3.587408148320037*^-15
    "Final Residual Error for Psi: "7.260858581048524*^-14
    "Final Residual Error for Psi: "1.3401314834324968*^-15
    "Final Residual Error for Psi: "1.3655526362454928*^-15
    "Final Residual Error for Psi: "1.4048346362073583*^-15
    "Final Residual Error for Psi: "4.7878367936959876*^-15
    "Final Residual Error for Psi: "3.566967549236896*^-15
    "Final Residual Error for Psi: "1.9047805867311585*^-15
    "Final Residual Error for Psi: "2.131975151975496*^-15
    >>> Print[stm1];
    {15.908635587619125, 788.4836843673284, 63966.1153877612}
    >>> Print[std1];
    {0.24787373186381415, 0.4464856780850617, 0.5791926324418597}
    >>> Print[stm2];
    {5.373495426968223, 102.74716094297389, 3650.989513270755}
    >>> Print[std2];
    {0.3213389051746825, 0.7054814215110067, 0.8391121876950028}
    >>> Print[stm3];
    {2.2551683685147204, 13.11393510238892, 124.72993277593807}
    >>> Print[std3];
    {0.45672177440186557, 0.8698892141420411, 0.9722166488850249}

    For Python/Numpy:

    >>> D0 = ml.matrix([[-5.49, 0., 1.15, 0.],[0., -2.29, 0., 0.],[0., 0.08, -1.32, 0.],[0.72, 1.17, 0.7, -7.07]])
    >>> D1 = ml.matrix([[0.25, 0.38, 0.64, 0.],[0., 0., 0., 1.09],[0., 1.24, 0., 0.],[0.37, 0., 0., 0.]])
    >>> D2 = ml.matrix([[0.3, 1.0, 0., 0.48],[0., 0.2, 0., 0.],[0., 0., 0., 0.],[0.61, 0., 0., 0.2]])
    >>> D3 = ml.matrix([[0., 0.98, 0., 0.31],[0., 0., 1.0, 0.],[0., 0., 0., 0.],[1.1, 0.84, 0.33, 1.03]])
    >>> sigma3 = ml.matrix([[0.83333,0.11404,0.05263]])
    >>> S3 = ml.matrix([[-3., 0., 0.],[0.73077, -0.73077, 0.],[0., 0.5, -0.5]])
    >>> sigma2 = ml.matrix([[1.]])
    >>> S2 = ml.matrix([[-2.]])
    >>> sigma1 = ml.matrix([[0.25,0.75]])
    >>> S1 = ml.matrix([[-2.5, 2.5],[0., -10.]])
    >>> ncm1, ncd1, ncm2, ncd2, ncm3, ncd3 = MMAPPH1NPPR([D0, D1, D2, D3], [sigma1, sigma2, sigma3], [S1, S2, S3], "ncMoms", 3, "ncDistr", 500)
    Final Residual Error for G:  2.3314683517128287e-15
    Final Residual Error for G:  2.3314683517128287e-15
    Final Residual Error for G:  1.790234627208065e-15
    Final Residual Error for G:  3.6914915568786455e-15
    Final Residual Error for G:  1.4521717162097048e-13
    Final Residual Error for G:  5.974387651264124e-15
    Final Residual Error for G:  1.0630385460785874e-14
    Final Residual Error for G:  2.7131075164277263e-15
    Final Residual Error for G:  3.20949754572e-16
    Final Residual Error for R:  3.63464120372e-16
    Final Residual Error for G:  3.20949754572e-16
    Final Residual Error for R:  3.63464120372e-16
    >>> distrPoints = [1., 5., 10.]
    >>> stm1, std1, stm2, std2, stm3, std3 = MMAPPH1NPPR([D0, D1, D2, D3], [sigma1, sigma2, sigma3], [S1, S2, S3], "stMoms", 3, "stDistr", distrPoints)
    Final Residual Error for G:  2.3314683517128287e-15
    Final Residual Error for G:  2.3314683517128287e-15
    Final Residual Error for G:  1.790234627208065e-15
    Final Residual Error for G:  3.6914915568786455e-15
    Final Residual Error for G:  1.4521717162097048e-13
    Final Residual Error for G:  9.319276714258098e-16
    Final Residual Error for G:  7.064661355915547e-16
    Final Residual Error for G:  9.459663954936026e-16
    Final Residual Error for G:  1.0630385460785874e-14
    Final Residual Error for G:  3.387294073432367e-15
    Final Residual Error for G:  2.328351270466933e-15
    Final Residual Error for G:  1.0900568642169262e-15
    >>> print(stm1)
    [15.908635587617834, 788.48368436722012, 63966.115387749385]
    >>> print(std1)
    [ 0.24787  0.44649  0.57919]
    >>> print(stm2)
    [5.3734954269680575, 102.74716094296781, 3650.9895132704323]
    >>> print(std2)
    [ 0.32134  0.70548  0.83911]
    >>> print(stm3)
    [2.2551683685147075, 13.113935102388783, 124.72993277593639]
    >>> print(std3)
    [ 0.45672  0.86989  0.97222]

