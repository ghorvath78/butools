butools.queues.MMAPPH1FCFS
==========================

.. currentmodule:: butools.queues

.. np:function:: MMAPPH1FCFS

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Ret = MMAPPH1FCFS(D, sigma, S, ...)`
        * - Mathematica:
          - :code:`Ret = MMAPPH1FCFS[D, sigma, S, ...]`
        * - Python/Numpy:
          - :code:`Ret = MMAPPH1FCFS(D, sigma, S, ...)`

    Returns various performane measures of a MMAP[K]/PH[K]/1 
    first-come-first-serve queue, see [1]_.

    Parameters
    ----------
    D : list of matrices of shape (N,N), length (K+1)
        The D0...DK matrices of the arrival process.
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
        | "stDistrME"    | None               | The vector-matrix parameters of the    |
        |                |                    | matrix-exponentially distributed       |
        |                |                    | sojourn time distribution              |
        +----------------+--------------------+----------------------------------------+
        | "stDistrPH"    | None               | The vector-matrix parameters of the    |
        |                |                    | matrix-exponentially distributed       |
        |                |                    | sojourn time distribution, converted   |
        |                |                    | to a continuous PH representation      |
        +----------------+--------------------+----------------------------------------+
        | "prec"         | The precision      | Numerical precision used as a stopping |
        |                |                    | condition when solving the Riccati     |
        |                |                    | equation                               |
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
    .. [1] Qiming He, "Analysis of a continuous time 
           SM[K]/PH[K]/1/FCFS queue: Age process, sojourn times,
           and queue lengths", Journal of Systems Science and 
           Complexity, 25(1), pp 133-155, 2012.

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
    >>> [ncm1, ncd1, ncm2, ncd2, ncm3, ncd3] = MMAPPH1FCFS({D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, 'ncMoms', 3, 'ncDistr', 500);
    Final Residual Error for Psi:    4.0523e-15
    >>> distrPoints = [1., 5., 10.];
    >>> [stm1, std1, stm2, std2, stm3, std3] = MMAPPH1FCFS({D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, 'stMoms', 3, 'stDistr', distrPoints);
    Final Residual Error for Psi:    4.0523e-15
    >>> disp(stm1);
            5.847       84.706       1866.7
    >>> disp(std1);
          0.28789      0.60379      0.79933
    >>> disp(stm2);
           6.3613       91.529       2014.4
    >>> disp(std2);
          0.20213      0.57229       0.7835
    >>> disp(stm3);
           6.4108       92.984       2049.2
    >>> disp(std3);
          0.21755       0.5651      0.77972
    >>> [betap1, Bp1, betap2, Bp2, betap3, Bp3] = MMAPPH1FCFS({D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, 'stDistrPH');
    Final Residual Error for Psi:    4.0523e-15
    >>> [beta1, B1, beta2, B2, beta3, B3] = MMAPPH1FCFS({D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, 'stDistrME');
    Final Residual Error for Psi:    4.0523e-15
    >>> stdFromPH1 = CdfFromPH(betap1, Bp1, distrPoints);
    >>> disp(stdFromPH1);
          0.28789      0.60379      0.79933
    >>> stmFromME1 = MomentsFromME(beta1, B1, 3);
    >>> disp(stmFromME1);
            5.847       84.706       1866.7
    >>> stdFromPH2 = CdfFromPH(betap2, Bp2, distrPoints);
    >>> disp(stdFromPH2);
          0.20213      0.57229       0.7835
    >>> stmFromME2 = MomentsFromME(beta2, B2, 3);
    >>> disp(stmFromME2);
           6.3613       91.529       2014.4
    >>> stdFromPH3 = CdfFromPH(betap3, Bp3, distrPoints);
    >>> disp(stdFromPH3);
          0.21755       0.5651      0.77972
    >>> stmFromME3 = MomentsFromME(beta3, B3, 3);
    >>> disp(stmFromME3);
           6.4108       92.984       2049.2

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
    >>> {ncm1, ncd1, ncm2, ncd2, ncm3, ncd3} = MMAPPH1FCFS[{D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, "ncMoms", 3, "ncDistr", 500];
    "Final Residual Error for Psi: "7.105427357601002*^-15
    >>> distrPoints = {1., 5., 10.};
    >>> {stm1, std1, stm2, std2, stm3, std3} = MMAPPH1FCFS[{D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, "stMoms", 3, "stDistr", distrPoints];
    "Final Residual Error for Psi: "7.105427357601002*^-15
    >>> Print[stm1];
    {5.847006216024009, 84.70619836863656, 1866.7093911950303}
    >>> Print[std1];
    {0.28789222215503396, 0.6037889208096201, 0.7993288816416165}
    >>> Print[stm2];
    {6.3613219350713415, 91.52937282421111, 2014.4076848738787}
    >>> Print[std2];
    {0.20213123601937755, 0.5722851898998502, 0.7834982143099977}
    >>> Print[stm3];
    {6.410840761114221, 92.98444287168671, 2049.1524939597166}
    >>> Print[std3];
    {0.21755284454422186, 0.5650994743810849, 0.7797186180777291}
    >>> {betap1, Bp1, betap2, Bp2, betap3, Bp3} = MMAPPH1FCFS[{D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, "stDistrPH"];
    "Final Residual Error for Psi: "7.105427357601002*^-15
    >>> {beta1, B1, beta2, B2, beta3, B3} = MMAPPH1FCFS[{D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, "stDistrME"];
    "Final Residual Error for Psi: "7.105427357601002*^-15
    >>> stdFromPH1 = CdfFromPH[betap1, Bp1, distrPoints];
    >>> Print[stdFromPH1];
    {0.2878922221550334, 0.6037889208096203, 0.7993288816416169}
    >>> stmFromME1 = MomentsFromME[beta1, B1, 3];
    >>> Print[stmFromME1];
    {5.847006216024715, 84.70619836865984, 1866.7093911957213}
    >>> stdFromPH2 = CdfFromPH[betap2, Bp2, distrPoints];
    >>> Print[stdFromPH2];
    {0.2021312360193771, 0.5722851898998502, 0.7834982143099979}
    >>> stmFromME2 = MomentsFromME[beta2, B2, 3];
    >>> Print[stmFromME2];
    {6.361321935071883, 91.52937282422667, 2014.4076848742416}
    >>> stdFromPH3 = CdfFromPH[betap3, Bp3, distrPoints];
    >>> Print[stdFromPH3];
    {0.21755284454422208, 0.5650994743810853, 0.7797186180777295}
    >>> stmFromME3 = MomentsFromME[beta3, B3, 3];
    >>> Print[stmFromME3];
    {6.410840761113892, 92.98444287168762, 2049.152493959693}

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
    >>> ncm1, ncd1, ncm2, ncd2, ncm3, ncd3 = MMAPPH1FCFS([D0, D1, D2, D3], [sigma1, sigma2, sigma3], [S1, S2, S3], "ncMoms", 3, "ncDistr", 500)
    Final Residual Error for G:  4.121702978920894e-15
    >>> distrPoints = [1., 5., 10.]
    >>> stm1, std1, stm2, std2, stm3, std3 = MMAPPH1FCFS([D0, D1, D2, D3], [sigma1, sigma2, sigma3], [S1, S2, S3], "stMoms", 3, "stDistr", distrPoints)
    Final Residual Error for G:  4.121702978920894e-15
    >>> print(stm1)
    [5.8470062160235257, 84.706198368623802, 1866.7093911946247]
    >>> print(std1)
    [ 0.28789  0.60379  0.79933]
    >>> print(stm2)
    [6.3613219350708565, 91.529372824197793, 2014.4076848734521]
    >>> print(std2)
    [ 0.20213  0.57229  0.7835 ]
    >>> print(stm3)
    [6.4108407611137332, 92.984442871673338, 2049.1524939592864]
    >>> print(std3)
    [ 0.21755  0.5651   0.77972]
    >>> betap1, Bp1, betap2, Bp2, betap3, Bp3 = MMAPPH1FCFS([D0, D1, D2, D3], [sigma1, sigma2, sigma3], [S1, S2, S3], "stDistrPH")
    Final Residual Error for G:  4.121702978920894e-15
    >>> beta1, B1, beta2, B2, beta3, B3 = MMAPPH1FCFS([D0, D1, D2, D3], [sigma1, sigma2, sigma3], [S1, S2, S3], "stDistrME")
    Final Residual Error for G:  4.121702978920894e-15
    >>> stdFromPH1 = CdfFromPH(betap1, Bp1, distrPoints)
    >>> print(stdFromPH1)
    [ 0.28789  0.60379  0.79933]
    >>> stmFromME1 = MomentsFromME(beta1, B1, 3)
    >>> print(stmFromME1)
    [5.847006216026557, 84.7061983687164, 1866.7093911974439]
    >>> stdFromPH2 = CdfFromPH(betap2, Bp2, distrPoints)
    >>> print(stdFromPH2)
    [ 0.20213  0.57229  0.7835 ]
    >>> stmFromME2 = MomentsFromME(beta2, B2, 3)
    >>> print(stmFromME2)
    [6.3613219350690784, 91.52937282416201, 2014.4076848721443]
    >>> stdFromPH3 = CdfFromPH(betap3, Bp3, distrPoints)
    >>> print(stdFromPH3)
    [ 0.21755  0.5651   0.77972]
    >>> stmFromME3 = MomentsFromME(beta3, B3, 3)
    >>> print(stmFromME3)
    [6.4108407611132225, 92.984442871662935, 2049.1524939590845]

