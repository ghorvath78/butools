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

        +----------------+--------------------+--------------------------------------+
        | Parameter name | Input parameters   | Output                               |
        +================+====================+======================================+
        | "qlMoms"       | Number of moments  | The queue length moments             |
        +----------------+--------------------+--------------------------------------+
        | "qlDistr"      | Upper bound (int)  | The queue length distribution is     |
        |                |                    | returned up to this bound.           |
        +----------------+--------------------+--------------------------------------+
        | "stMoms"       | Number of moments  | The sojourn time moments             |
        +----------------+--------------------+--------------------------------------+
        | "stDistr"      | A vector of points | The sojourn time distribution at the |
        |                |                    | requested points (cummulative, cdf)  |
        +----------------+--------------------+--------------------------------------+
        | "prec"         | The precision      | Numerical precision to check if the  |
        |                |                    | input is valid and it is also used   |
        |                |                    | as a stopping condition when solving |
        |                |                    | the matrix-quadratic equation        |
        +----------------+--------------------+--------------------------------------+
        | "erlMaxOrder"  | Integer number     | The maximal Erlang order used in the |
        |                |                    | erlangization procedure. The default |
        |                |                    | value is 200.                        |
        +----------------+--------------------+--------------------------------------+
        | "classes"      | Vector of integers | Only the performance measures        |
        |                |                    | belonging to these classes are       |
        |                |                    | returned. If not given, all classes  |
        |                |                    | are analyzed.                        |
        +----------------+--------------------+--------------------------------------+

        (The queue length related quantities include the customer 
        in the server, and the sojourn time related quantities 
        include the service times as well)

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
           Research, 2015, to appear.
           doi:10.1016/j.ejor.2015.03.004

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
    >>> [qlm, qld] = MMAPPH1PRPR({D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, 'qlMoms', 3, 'qlDistr', 500);
    Final Residual Error for Psi:     4.774e-15
    Final Residual Error for Psi:    6.7946e-14
    Final Residual Error for Psi:    5.9258e-15
    Final Residual Error for Psi:    3.8719e-15
    Final Residual Error for Psi:    7.1054e-15
    Final Residual Error for Psi:     2.498e-15
    Final Residual Error for Psi:    2.7478e-15
    >>> distrPoints = [1., 5., 10.];
    >>> [stm, std] = MMAPPH1PRPR({D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, 'stMoms', 3, 'stDistr', distrPoints);
    Final Residual Error for Psi:     4.774e-15
    Final Residual Error for Psi:    6.7946e-14
    Final Residual Error for Psi:     3.804e-15
    Final Residual Error for Psi:    3.6332e-15
    Final Residual Error for Psi:    3.2784e-15
    Final Residual Error for Psi:    3.8719e-15
    Final Residual Error for Psi:    7.1054e-15
    Final Residual Error for Psi:    5.3424e-16
    Final Residual Error for Psi:    8.5008e-16
    Final Residual Error for Psi:    9.6955e-16
    Final Residual Error for Psi:    2.7478e-15
    >>> disp(stm);
           16.306       5.5053       1.9703
           806.63       107.78       11.741
            65397         3848       112.85
    >>> disp(std);
          0.22164      0.34346      0.55657
          0.43134      0.69324      0.88197
          0.56865      0.83131      0.97479

    For Mathematica:

    
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
    >>> qlm, qld = MMAPPH1PRPR([D0, D1, D2, D3], [sigma1, sigma2, sigma3], [S1, S2, S3], "qlMoms", 3, "qlDistr", 500)
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
    >>> stm, std = MMAPPH1PRPR([D0, D1, D2, D3], [sigma1, sigma2, sigma3], [S1, S2, S3], "stMoms", 3, "stDistr", distrPoints)
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
    >>> print(stm)
    [[  1.63060e+01   5.50526e+00   1.97033e+00]
     [  8.06630e+02   1.07777e+02   1.17412e+01]
     [  6.53972e+04   3.84804e+03   1.12849e+02]]
    >>> print(std)
    [[ 0.22164  0.34346  0.55657]
     [ 0.43134  0.69324  0.88197]
     [ 0.56865  0.83131  0.97479]]

