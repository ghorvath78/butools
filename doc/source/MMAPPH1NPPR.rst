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
    For MATLAB:

    >>> D0=[-5.49, 0, 1.15, 0; 0, -2.29, 0, 0; 0, 0.08, -1.32, 0; 0.72, 1.17, 0.7, -7.07];
    >>> D1=[0.25, 0.38, 0.64, 0; 0, 0, 0, 1.09; 0, 1.24, 0, 0; 0.37, 0, 0, 0];
    >>> D2=[0.3, 1.0, 0, 0.48; 0, 0.2, 0, 0; 0, 0, 0, 0; 0.61, 0, 0, 0.2];
    >>> D3=[0, 0.98, 0, 0.31; 0, 0, 1.0, 0; 0, 0, 0, 0; 1.1, 0.84, 0.33, 1.03];
    >>> sigma1 = [1/4, 1-1/4];
    >>> S1 = [-2.5, 2.5; 0, -10];
    >>> sigma2 = [1/1.7, 1-1/1.7];
    >>> S2 = [-2/0.85, 2/0.85; 0, -4];
    >>> sigma3 = [1/6, 5/6];
    >>> S3 = [-0.5, 0.5; 0, -3];
    >>> [stm, std] = MMAPPH1NPPR({D0,D1,D2,D3}, {sigma1,sigma2,sigma3},{S1,S2,S3}, 'stMoms', 3, 'stDistr', [1,5,10]);
    >>> std
          0.24813       0.3202      0.45966
          0.44751      0.70733      0.87042
          0.58038      0.84018      0.97233
    >>> stm
           15.831       5.3482       2.2458
           781.09       101.93       13.054
            63091       3617.3       124.21
    >>> meanQlFromLittlesLaw = [stm(1,1)/MarginalMomentsFromMAP(D0+D2+D3,D1,1), stm(1,2)/MarginalMomentsFromMAP(D0+D1+D3,D2,1), stm(1,3)/MarginalMomentsFromMAP(D0+D1+D2,D3,1)]
           17.492       1.2717       1.7072
    >>> [qlm, qld] = MMAPPH1NPPR({D0,D1,D2,D3}, {sigma1,sigma2,sigma3},{S1,S2,S3}, 'qlMoms', 3, 'qlDistr', 50);
    >>> qlm
           17.931       1.3035       1.4978
             1022       7.9555       8.1379
            95668       82.606       68.996
    >>> plot(qld)
  
    For Python/Numpy:
    
    >>> D0=ml.matrix([[-5.49, 0, 1.15, 0], [ 0, -2.29, 0, 0], [ 0, 0.08, -1.32, 0], [ 0.72, 1.17, 0.7, -7.07]])
    >>> D1=ml.matrix([[0.25, 0.38, 0.64, 0], [ 0, 0, 0, 1.09], [ 0, 1.24, 0, 0], [ 0.37, 0, 0, 0]])
    >>> D2=ml.matrix([[0.3, 1.0, 0, 0.48], [ 0, 0.2, 0, 0], [ 0, 0, 0, 0], [ 0.61, 0, 0, 0.2]])
    >>> D3=ml.matrix([[0, 0.98, 0, 0.31], [ 0, 0, 1.0, 0], [ 0, 0, 0, 0], [ 1.1, 0.84, 0.33, 1.03]])
    >>> sigma3 = ml.matrix([[1/6, 5/6]])
    >>> S3 = ml.matrix([[-0.5, 0.5], [ 0, -3]])
    >>> sigma2 = ml.matrix([[1/1.7, 1-1/1.7]])
    >>> S2 = ml.matrix([[-2/0.85, 2/0.85], [ 0, -4]])
    >>> sigma1 = ml.matrix([[1/4, 1-1/4]])
    >>> S1 = ml.matrix([[-2.5, 2.5], [ 0, -10]])
    >>> stm, std = MMAPPH1NPPR((D0,D1,D2,D3), (sigma1,sigma2,sigma3),(S1,S2,S3), 'stMoms', 3, 'stDistr', [1,5,10]) 
    >>> print(std)
    [[ 0.24813421  0.32020161  0.459657  ]
     [ 0.44751045  0.70733126  0.87041572]
     [ 0.58037807  0.84018473  0.97232579]]
    >>> print(stm)
    [[  1.58305010e+01   5.34816990e+00   2.24576568e+00]
     [  7.81090522e+02   1.01931494e+02   1.30540555e+01]
     [  6.30907058e+04   3.61733345e+03   1.24208816e+02]]
    >>> meanQlFromLittlesLaw = np.array([stm[0,0]/MarginalMomentsFromMAP(D0+D2+D3,D1,1)[0], stm[0,1]/MarginalMomentsFromMAP(D0+D1+D3,D2,1)[0], stm[0,2]/MarginalMomentsFromMAP(D0+D1+D2,D3,1)[0]])
    [ 17.49185168   1.2716964    1.70720665]
    >>> print(meanQlFromLittlesLaw)
    >>> qlm, qld = MMAPPH1NPPR((D0,D1,D2,D3), (sigma1,sigma2,sigma3),(S1,S2,S3), 'qlMoms', 3, 'qlDistr', 50)
    >>> print(qlm)
    [[  1.74918517e+01   1.27169640e+00   1.70720665e+00]
     [  9.98469939e+02   7.60714872e+00   9.08782742e+00]
     [  9.35269403e+04   7.84157193e+01   7.62178217e+01]]
    >>> plt.plot(qld)

