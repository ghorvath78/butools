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
    >>> [stm, std] = MMAPPH1PRPR({D0,D1,D2,D3}, {sigma1,sigma2,sigma3},{S1,S2,S3}, 'stMoms', 3, 'stDistr', [1,5,10]);
    >>> std
          0.22167      0.34046      0.55657
          0.43232      0.69522      0.88197
          0.56982      0.83254      0.97479
    >>> stm
           16.228        5.482       1.9703
           799.15       106.86       11.741
            64508       3809.4       112.85
    >>> meanQlFromLittlesLaw = [stm(1,1)/MarginalMomentsFromMAP(D0+D2+D3,D1,1), stm(1,2)/MarginalMomentsFromMAP(D0+D1+D3,D2,1), stm(1,3)/MarginalMomentsFromMAP(D0+D1+D2,D3,1)]
           17.931       1.3035       1.4978
    >>> [qlm, qld] = MMAPPH1PRPR({D0,D1,D2,D3}, {sigma1,sigma2,sigma3},{S1,S2,S3}, 'qlMoms', 3, 'qlDistr', 50);
    >>> qlm
           17.492       1.2717       1.7072
           998.47       7.6071       9.0878
            93527       78.416       76.218
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
    >>> stm, std = MMAPPH1PRPR((D0,D1,D2,D3), (sigma1,sigma2,sigma3),(S1,S2,S3), 'stMoms', 3, 'stDistr', [1,5,10]) 
    >>> print(std)
    [[ 0.22166892  0.34045801  0.55657091]
     [ 0.43231964  0.69522131  0.88196947]
     [ 0.56982261  0.83253605  0.97479177]]
    >>> print(stm)
    [[  1.62280212e+01   5.48201294e+00   1.97033167e+00]
     [  7.99146185e+02   1.06860596e+02   1.17412894e+01]
     [  6.45081685e+04   3.80943440e+03   1.12850305e+02]]
    >>> meanQlFromLittlesLaw = np.array([stm[0,0]/MarginalMomentsFromMAP(D0+D2+D3,D1,1)[0], stm[0,1]/MarginalMomentsFromMAP(D0+D1+D3,D2,1)[0], stm[0,2]/MarginalMomentsFromMAP(D0+D1+D2,D3,1)[0]])
    [ 17.9310901    1.30352182   1.49782471]
    >>> print(meanQlFromLittlesLaw)
    >>> qlm, qld = MMAPPH1PRPR((D0,D1,D2,D3), (sigma1,sigma2,sigma3),(S1,S2,S3), 'qlMoms', 3, 'qlDistr', 50)
    >>> print(qlm)
    [[  1.79310901e+01   1.30352182e+00   1.49782471e+00]
     [  1.02199118e+03   7.95547352e+00   8.13794999e+00]
     [  9.56678067e+04   8.26061428e+01   6.89963754e+01]]
    >>> plt.plot(qld)

