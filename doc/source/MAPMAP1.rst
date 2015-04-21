butools.queues.MAPMAP1
======================

.. currentmodule:: butools.queues

.. np:function:: MAPMAP1

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Ret = MAPMAP1(D0, D1, S0, S1, ...)`
        * - Mathematica:
          - :code:`Ret = MAPMAP1[D0, D1, S0, S1, ...]`
        * - Python/Numpy:
          - :code:`Ret = MAPMAP1(D0, D1, S0, S1, ...)`

    Returns various performane measures of a MAP/MAP/1 queue.

    In a MAP/MAP/1 queue both the arrival and the service
    processes are characterized by Markovian arrival 
    processes.

    Parameters
    ----------
    D0 : matrix, shape(N,N)
        The transitions of the arrival MAP not accompanied by
        job arrivals
    D1 : matrix, shape(N,N)
        The transitions of the arrival MAP accompanied by
        job arrivals
    S0 : matrix, shape(N,N)
        The transitions of the service MAP not accompanied by
        job service
    S1 : matrix, shape(N,N)
        The transitions of the service MAP accompanied by
        job service
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
        | "qlDistr"      | A vector of points | The queue length distribution at     |
        |                |                    | the requested points                 |
        +----------------+--------------------+--------------------------------------+
        | "qlDistrMG"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-geometrically distributed     |
        |                |                    | queue length distribution            |
        +----------------+--------------------+--------------------------------------+
        | "qlDistrDPH"   | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-geometrically distributed     |
        |                |                    | queue length distribution, converted |
        |                |                    | to a discrete PH representation      |
        +----------------+--------------------+--------------------------------------+
        | "stMoms"       | Number of moments  | The sojourn time moments             |
        +----------------+--------------------+--------------------------------------+
        | "stDistr"      | A vector of points | The sojourn time distribution at the |
        |                |                    | requested points (cummulative, cdf)  |
        +----------------+--------------------+--------------------------------------+
        | "stDistrME"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-exponentially distributed     |
        |                |                    | sojourn time distribution            |
        +----------------+--------------------+--------------------------------------+
        | "stDistrPH"    | None               | The vector-matrix parameters of the  |
        |                |                    | matrix-exponentially distributed     |
        |                |                    | sojourn time distribution, converted |
        |                |                    | to a continuous PH representation    |
        +----------------+--------------------+--------------------------------------+
        | "prec"         | The precision      | Numerical precision to check if the  |
        |                |                    | input is valid and it is also used   |
        |                |                    | as a stopping condition when solving |
        |                |                    | the matrix-quadratic equation        |
        +----------------+--------------------+--------------------------------------+

        (The queue length related quantities include the customer 
        in the server, and the sojourn time related quantities 
        include the service times as well)

    Returns
    -------
    Ret : list of the performance measures
        Each entry of the list corresponds to a performance 
        measure requested. If there is just a single item, 
        then it is not put into a list.
   
    Notes
    -----
    "qlDistrMG" and "stDistrME" behave much better numerically than 
    "qlDistrDPH" and "stDistrPH".

    Examples
    ========    
    For MATLAB:

    >>> D0 = [-8, 1, 2; 0, -6, 4; 3 0 -3]
    >>> D1 = [4, 1, 0; 0, 2, 0; 0, 0, 0]
    >>> S0=[-10, 4; 0, -7]
    >>> S1=[5, 1; 4, 3]
    >>> [qld,qlm] = MAPMAP1(D0, D1, S0, S1, 'qlDistr', (0:10), 'qlMoms', 5);
    >>> qld
          0.67697      0.18891      0.07951     0.032563     0.013182    0.0053087    0.0021328   0.00085583   0.00034322    0.0001376   5.5158e-05
    >>> qlm
          0.54864        1.306        4.357       19.193       105.39
    >>> [alphap,Ap] = MAPMAP1(D0, D1, S0, S1, 'qlDistrDPH');
    >>> alphap
         0.067133     0.075328       0.0407     0.047054     0.035736     0.057079
    >>> Ap
           0.2813     0.046292    0.0076126    0.0034456            0            0
          0.10552      0.33506    0.0083548      0.01223            0            0
          0.17544     0.045811      0.13666     0.010563            0            0
         0.079843      0.20353     0.042265      0.16552            0            0
          0.20818     0.091826      0.05939     0.013516            0            0
           0.1325      0.20874     0.034792     0.065764            0            0
    >>> PmfFromDPH(alphap,Ap,(0:10))'                
          0.67697      0.18891      0.07951     0.032563     0.013182    0.0053087    0.0021328   0.00085583   0.00034322    0.0001376   5.5158e-05
    >>> MomentsFromDPH(alphap,Ap,5)
          0.54864        1.306        4.357       19.193       105.39
    >>> [std, stm] = MAPMAP1(D0, D1, S0, S1, 'stDistr', (0:0.1:1), 'stMoms', 5);
    >>> std
       1.1102e-16       0.3122      0.53197      0.68345      0.78667      0.85655      0.90367      0.93538      0.95667      0.97097      0.98055
    >>> stm
          0.25908      0.13145      0.09911     0.099178      0.12376
    >>> [beta, B] = MAPMAP1(D0, D1, S0, S1, 'stDistrME');
    >>> beta
          0.40071      0.26596      0.19675      0.13659            0            0
    >>> B
          -8.1822        4.607      0.61259      0.21807            0            0
           1.4951      -5.9854      0.51299      0.33032            0            0
          0.24927      0.11789      -9.0589       4.3127            0            0
          0.20579       0.1273      0.75919      -6.4625            0            0
           0.5575      0.23374      0.18872     0.080814          -10            4
          0.44853      0.30439       0.1539     0.099095            0           -7
    >>> CdfFromME(beta,B,(0:0.1:1))'   
       1.1102e-16       0.3122      0.53197      0.68345      0.78667      0.85655      0.90367      0.93538      0.95667      0.97097      0.98055
    >>> MomentsFromME(beta,B,5)
          0.25908      0.13145      0.09911     0.099178      0.12376
  

