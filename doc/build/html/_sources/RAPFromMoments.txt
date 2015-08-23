butools.map.RAPFromMoments
==========================

.. currentmodule:: butools.map

.. np:function:: RAPFromMoments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[H0, H1] = RAPFromMoments(moms, Nm)`
        * - Mathematica:
          - :code:`{H0, H1} = RAPFromMoments[moms, Nm]`
        * - Python/Numpy:
          - :code:`H0, H1 = RAPFromMoments(moms, Nm)`

    Creates a rational arrival process that has the same 
    marginal and lag-1 joint moments as given (see [1]_).

    Parameters
    ----------
    moms : vector of doubles
        The list of marginal moments. To obtain a rational 
        process of order M, 2*M-1 marginal moments are 
        required.
    Nm : matrix, shape (M,M)
        The matrix of lag-1 joint moments. 
    
    Returns
    -------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational process
    
    Notes
    -----
    There is no guarantee that the returned matrices define
    a valid stochastic process. The joint densities may be
    negative.

    References
    ----------
    .. [1] G Horvath, M Telek, "A minimal representation of 
           Markov arrival processes and a moments matching 
           method," Performance Evaluation 64:(9-12) pp. 
           1153-1168. (2007)       

    Examples
    ========
    For Matlab:

    >>> G0 = [-6.2, 2., 0.; 2., -9., 1.; 1., 0, -3.];
    >>> G1 = [2.2, -2., 4.; 2., 2., 2.; 1., 0, 1.];
    >>> moms = MarginalMomentsFromRAP(G0, G1, 5);
    >>> disp(moms);
          0.36585      0.25535      0.26507      0.36691      0.63573
    >>> Nm = LagkJointMomentsFromRAP(G0, G1, 2, 1);
    >>> disp(Nm);
                1      0.36585      0.25535
          0.36585      0.12866     0.088334
          0.25535     0.088802      0.06067
    >>> [H0, H1] = RAPFromMoments(moms, Nm);
    >>> disp(H0);
          -12.949        36.78      -24.817
          -1.1102      -2.5113      0.91705
         -0.71205      0.68912      -2.7393
    >>> disp(H1);
           9.2672      -99.958       91.678
           1.1693      -2.1771       3.7123
          0.65292       3.9994      -1.8901
    >>> rmoms = MarginalMomentsFromRAP(H0, H1, 5);
    >>> disp(rmoms);
          0.36585      0.25535      0.26507      0.36691      0.63573
    >>> rNm = LagkJointMomentsFromRAP(H0, H1, 2, 1);
    >>> disp(rNm);
                1      0.36585      0.25535
          0.36585      0.12866     0.088334
          0.25535     0.088802      0.06067
    >>> G0 = [-5., 0, 1., 1.; 1., -8., 1., 0; 1., 0, -4., 1.; 1., 2., 3., -9.];
    >>> G1 = [0, 1., 0, 2.; 2., 1., 3., 0; 0, 0, 1., 1.; 1., 1., 0, 1.];
    >>> moms = MarginalMomentsFromRAP(G0, G1, 7);
    >>> disp(moms);
      Columns 1 through 6
          0.34247      0.25054      0.28271      0.42984      0.81999       1.8795
      Column 7
            5.028
    >>> Nm = LagkJointMomentsFromRAP(G0, G1, 3, 1);
    >>> disp(Nm);
                1      0.34247      0.25054      0.28271
          0.34247       0.1173     0.085789     0.096807
          0.25054       0.0857     0.062633      0.07066
          0.28271     0.096627     0.070589     0.079623
    >>> [H0, H1] = RAPFromMoments(moms, Nm);
    >>> disp(H0);
          -6.7126       32.989      -108.71       77.983
          -0.8704      -8.3405       25.268      -19.013
         -0.65982       3.0739      -16.543       11.227
         -0.65977       3.0766      -10.915       5.5959
    >>> disp(H1);
           1.6406       4.4202        26351       -26353
          0.61635      0.87632      -1431.5       1432.9
          0.78683      0.65689       716.34      -714.88
          0.78683      0.65679       717.32      -715.86
    >>> BuToolsCheckPrecision = 10.^-8;
    >>> rmoms = MarginalMomentsFromRAP(H0, H1, 7);
    >>> disp(rmoms);
      Columns 1 through 6
          0.34247      0.25054      0.28271      0.42984      0.81999       1.8795
      Column 7
            5.028
    >>> rNm = LagkJointMomentsFromRAP(H0, H1, 3, 1);
    >>> disp(rNm);
                1      0.34247      0.25054      0.28271
          0.34247       0.1173     0.085789     0.096807
          0.25054       0.0857     0.062633      0.07066
          0.28271     0.096627     0.070589     0.079623

    For Mathematica:

    
    For Python/Numpy:

    >>> G0 = ml.matrix([[-6.2, 2., 0.],[2., -9., 1.],[1., 0, -3.]])
    >>> G1 = ml.matrix([[2.2, -2., 4.],[2., 2., 2.],[1., 0, 1.]])
    >>> moms = MarginalMomentsFromRAP(G0, G1, 5)
    >>> print(moms)
    [0.36585365853658536, 0.25535027188603043, 0.26507255497329191, 0.36691170692675046, 0.635727559166956]
    >>> Nm = LagkJointMomentsFromRAP(G0, G1, 2, 1)
    >>> print(Nm)
    [[ 1.       0.36585  0.25535]
     [ 0.36585  0.12866  0.08833]
     [ 0.25535  0.0888   0.06067]]
    >>> H0, H1 = RAPFromMoments(moms, Nm)
    >>> print(H0)
    [[-12.94939  36.77983 -24.81743]
     [ -1.11017  -2.51134   0.91705]
     [ -0.71205   0.68912  -2.73927]]
    >>> print(H1)
    [[  9.26723 -99.95808  91.67785]
     [  1.16931  -2.17715   3.7123 ]
     [  0.65292   3.99937  -1.89008]]
    >>> rmoms = MarginalMomentsFromRAP(H0, H1, 5)
    >>> print(rmoms)
    [0.36585365853658869, 0.25535027188603343, 0.26507255497329529, 0.36691170692675534, 0.635727559166965]
    >>> rNm = LagkJointMomentsFromRAP(H0, H1, 2, 1)
    >>> print(rNm)
    [[ 1.       0.36585  0.25535]
     [ 0.36585  0.12866  0.08833]
     [ 0.25535  0.0888   0.06067]]
    >>> G0 = ml.matrix([[-5., 0, 1., 1.],[1., -8., 1., 0],[1., 0, -4., 1.],[1., 2., 3., -9.]])
    >>> G1 = ml.matrix([[0, 1., 0, 2.],[2., 1., 3., 0],[0, 0, 1., 1.],[1., 1., 0, 1.]])
    >>> moms = MarginalMomentsFromRAP(G0, G1, 7)
    >>> print(moms)
    [0.34246575342465752, 0.25053639214391815, 0.28270969431684256, 0.42984404959582057, 0.81998555487921787, 1.8794706921737607, 5.0280196843561171]
    >>> Nm = LagkJointMomentsFromRAP(G0, G1, 3, 1)
    >>> print(Nm)
    [[ 1.       0.34247  0.25054  0.28271]
     [ 0.34247  0.1173   0.08579  0.09681]
     [ 0.25054  0.0857   0.06263  0.07066]
     [ 0.28271  0.09663  0.07059  0.07962]]
    >>> H0, H1 = RAPFromMoments(moms, Nm)
    >>> print(H0)
    [[  -6.7126    32.98913 -108.70525   77.98325]
     [  -0.8704    -8.34051   25.2677   -19.01258]
     [  -0.65982    3.07388  -16.54281   11.22665]
     [  -0.65977    3.07663  -10.91488    5.59593]]
    >>> print(H1)
    [[  1.64061e+00   4.42021e+00   2.63514e+04  -2.63530e+04]
     [  6.16347e-01   8.76325e-01  -1.43147e+03   1.43293e+03]
     [  7.86827e-01   6.56888e-01   7.16343e+02  -7.14884e+02]
     [  7.86826e-01   6.56788e-01   7.17318e+02  -7.15860e+02]]
    >>> butools.checkPrecision = 10.**-8
    >>> rmoms = MarginalMomentsFromRAP(H0, H1, 7)
    >>> print(rmoms)
    [0.34246575341348912, 0.25053639213218126, 0.2827096943020489, 0.42984404957246941, 0.81998555483408853, 1.8794706920698447, 5.0280196840776616]
    >>> rNm = LagkJointMomentsFromRAP(H0, H1, 3, 1)
    >>> print(rNm)
    [[ 1.       0.34247  0.25054  0.28271]
     [ 0.34247  0.1173   0.08579  0.09681]
     [ 0.25054  0.0857   0.06263  0.07066]
     [ 0.28271  0.09663  0.07059  0.07962]]

