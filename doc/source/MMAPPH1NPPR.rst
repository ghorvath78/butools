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

    Returns various performane measures of a MMAP[K]/PH[K]/1 
    non-preemptive priority queue, see [1]_.

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
    >>> [qlm, qld] = MMAPPH1NPPR({D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, 'qlMoms', 3, 'qlDistr', 500);
    Final Residual Error for Psi:     4.774e-15
    Final Residual Error for Psi:     4.774e-15
    Final Residual Error for Psi:    3.0531e-15
    Final Residual Error for Psi:    3.5458e-15
    Final Residual Error for Psi:    7.7049e-14
    Final Residual Error for Psi:    2.5396e-15
    Final Residual Error for Psi:    7.1887e-15
    Final Residual Error for Psi:    4.5866e-15
    >>> distrPoints = [1., 5., 10.];
    >>> [stm, std] = MMAPPH1NPPR({D0, D1, D2, D3}, {sigma1, sigma2, sigma3}, {S1, S2, S3}, 'stMoms', 3, 'stDistr', distrPoints);
    Final Residual Error for Psi:     4.774e-15
    Final Residual Error for Psi:     4.774e-15
    Final Residual Error for Psi:    3.0531e-15
    Final Residual Error for Psi:    3.5458e-15
    Final Residual Error for Psi:    7.7049e-14
    Final Residual Error for Psi:     9.243e-16
    Final Residual Error for Psi:    1.1933e-15
    Final Residual Error for Psi:    1.3175e-15
    Final Residual Error for Psi:    7.1887e-15
    Final Residual Error for Psi:    3.7603e-15
    Final Residual Error for Psi:     2.319e-15
    Final Residual Error for Psi:    2.0864e-15
    >>> disp(stm);
           15.909       5.3735       2.2552
           788.48       102.75       13.114
            63966         3651       124.73
    >>> disp(std);
          0.24787      0.32134      0.45672
          0.44649      0.70548      0.86989
          0.57919      0.83911      0.97222

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
    >>> qlm, qld = MMAPPH1NPPR([D0, D1, D2, D3], [sigma1, sigma2, sigma3], [S1, S2, S3], "qlMoms", 3, "qlDistr", 500)
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
    >>> stm, std = MMAPPH1NPPR([D0, D1, D2, D3], [sigma1, sigma2, sigma3], [S1, S2, S3], "stMoms", 3, "stDistr", distrPoints)
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
    >>> print(stm)
    [[  1.59086e+01   5.37350e+00   2.25517e+00]
     [  7.88484e+02   1.02747e+02   1.31139e+01]
     [  6.39661e+04   3.65099e+03   1.24730e+02]]
    >>> print(std)
    [[ 0.24787  0.32134  0.45672]
     [ 0.44649  0.70548  0.86989]
     [ 0.57919  0.83911  0.97222]]

