butools.queues.MMAPPH1PRPR
==========================

.. currentmodule:: butools.queues

.. np:function:: MMAPPH1PRPR

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Ret = MMAPPH1PRPR(D, sigma, S, ...)`
        * - Mathematica:
          - :code:`Ret = MMAPPH1PRPR[D, sigma, S, ...]`
        * - Python/Numpy:
          - :code:`Ret = MMAPPH1PRPR(D, sigma, S, ...)`

    Returns various performane measures of a MMAP[K]/PH[K]/1 
    preemptive resume priority queue, see [1]_.

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
    >>> [ncm1, ncd1, ncm2, ncd2, ncm3, ncd3] = MMAPPH1PRPR({D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, 'ncMoms', 3, 'ncDistr', 500);
    Final Residual Error for Psi:     4.774e-15
    Final Residual Error for Psi:    6.8168e-14
    Final Residual Error for Psi:    5.9258e-15
    Final Residual Error for Psi:    3.8719e-15
    Final Residual Error for Psi:    6.5781e-15
    Final Residual Error for Psi:    2.5535e-15
    Final Residual Error for Psi:    2.7478e-15
    >>> distrPoints = [1., 5., 10.];
    >>> [stm1, std1, stm2, std2, stm3, std3] = MMAPPH1PRPR({D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, 'stMoms', 3, 'stDistr', distrPoints);
    Final Residual Error for Psi:     4.774e-15
    Final Residual Error for Psi:    6.8168e-14
    Final Residual Error for Psi:     3.804e-15
    Final Residual Error for Psi:    3.6351e-15
    Final Residual Error for Psi:    3.2784e-15
    Final Residual Error for Psi:    3.8719e-15
    Final Residual Error for Psi:    6.5781e-15
    Final Residual Error for Psi:    5.3424e-16
    Final Residual Error for Psi:    8.5008e-16
    Final Residual Error for Psi:    9.6955e-16
    Final Residual Error for Psi:    2.7478e-15
    >>> disp(stm1);
           16.306       806.63        65397
    >>> disp(std1);
          0.22164      0.43134      0.56865
    >>> disp(stm2);
           5.5053       107.78         3848
    >>> disp(std2);
          0.34346      0.69324      0.83131
    >>> disp(stm3);
           1.9703       11.741       112.85
    >>> disp(std3);
          0.55657      0.88197      0.97479

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
    >>> {ncm1, ncd1, ncm2, ncd2, ncm3, ncd3} = MMAPPH1PRPR[{D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, "ncMoms", 3, "ncDistr", 500];
    "Final Residual Error for Psi: "6.217248937900877*^-15
    "Final Residual Error for Psi: "6.864647739135421*^-14
    "Final Residual Error for Psi: "6.013418929473602*^-15
    "Final Residual Error for Psi: "3.8441472227646045*^-15
    "Final Residual Error for Psi: "7.022160630754115*^-15
    "Final Residual Error for Psi: "2.7200464103316335*^-15
    "Final Residual Error for Psi: "3.0600522116230877*^-15
    "Final Residual Error for G: "2.6259376617598917*^-16
    "Final Residual Error for R: "3.885780586188048*^-16
    "Final Residual Error for G: "2.6259376617598917*^-16
    "Final Residual Error for R: "3.885780586188048*^-16
    >>> distrPoints = {1., 5., 10.};
    >>> {stm1, std1, stm2, std2, stm3, std3} = MMAPPH1PRPR[{D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, "stMoms", 3, "stDistr", distrPoints];
    "Final Residual Error for Psi: "6.217248937900877*^-15
    "Final Residual Error for Psi: "6.864647739135421*^-14
    "Final Residual Error for Psi: "4.468647674116255*^-15
    "Final Residual Error for Psi: "4.649058915617843*^-15
    "Final Residual Error for Psi: "3.58046925441613*^-15
    "Final Residual Error for Psi: "3.8441472227646045*^-15
    "Final Residual Error for Psi: "7.022160630754115*^-15
    "Final Residual Error for Psi: "8.94618421124267*^-16
    "Final Residual Error for Psi: "1.818423883692688*^-15
    "Final Residual Error for Psi: "1.7907292944492614*^-15
    "Final Residual Error for Psi: "3.0600522116230877*^-15
    >>> Print[stm1];
    {16.30601164079545, 806.6298525859536, 65397.23294536998}
    >>> Print[std1];
    {0.2216364120343151, 0.4313440180099754, 0.5686546905139109}
    >>> Print[stm2];
    {5.505255217665003, 107.77744442363706, 3848.0354975573623}
    >>> Print[std2];
    {0.3434570011135701, 0.6932381350541303, 0.8313127970800466}
    >>> Print[stm3];
    {1.970327637987947, 11.741190396361796, 112.84868513581279}
    >>> Print[std3];
    {0.5565694682318661, 0.8819702634824608, 0.9747921567789963}

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
    >>> ncm1, ncd1, ncm2, ncd2, ncm3, ncd3 = MMAPPH1PRPR([D0, D1, D2, D3], [sigma1, sigma2, sigma3], [S1, S2, S3], "ncMoms", 3, "ncDistr", 500)
    Final Residual Error for G:  2.3314683517128287e-15
    Final Residual Error for G:  1.4521717162097048e-13
    Final Residual Error for G:  9.647230930776018e-15
    Final Residual Error for G:  1.5404344466674047e-15
    Final Residual Error for G:  6.7307270867900115e-15
    Final Residual Error for G:  1.7069679003611782e-15
    Final Residual Error for G:  1.339206523454095e-15
    Final Residual Error for G:  2.69315819645e-16
    Final Residual Error for R:  3.26128013484e-16
    Final Residual Error for G:  2.69315819645e-16
    Final Residual Error for R:  3.26128013484e-16
    >>> distrPoints = [1., 5., 10.]
    >>> stm1, std1, stm2, std2, stm3, std3 = MMAPPH1PRPR([D0, D1, D2, D3], [sigma1, sigma2, sigma3], [S1, S2, S3], "stMoms", 3, "stDistr", distrPoints)
    Final Residual Error for G:  2.3314683517128287e-15
    Final Residual Error for G:  1.4521717162097048e-13
    Final Residual Error for G:  7.88934883286745e-15
    Final Residual Error for G:  4.212775961409676e-15
    Final Residual Error for G:  3.1259717037102064e-15
    Final Residual Error for G:  1.5404344466674047e-15
    Final Residual Error for G:  6.7307270867900115e-15
    Final Residual Error for G:  7.642364031024785e-16
    Final Residual Error for G:  7.629124111965813e-16
    Final Residual Error for G:  4.657732532997727e-16
    Final Residual Error for G:  1.339206523454095e-15
    >>> print(stm1)
    [16.306011640793987, 806.62985258583274, 65397.232945356853]
    >>> print(std1)
    [ 0.22164  0.43134  0.56865]
    >>> print(stm2)
    [5.505255217664927, 107.77744442363372, 3848.0354975571686]
    >>> print(std2)
    [ 0.34346  0.69324  0.83131]
    >>> print(stm3)
    [1.9703276379879364, 11.741190396361684, 112.84868513581137]
    >>> print(std3)
    [ 0.55657  0.88197  0.97479]

