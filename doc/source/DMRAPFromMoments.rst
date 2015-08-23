butools.dmap.DMRAPFromMoments
=============================

.. currentmodule:: butools.dmap

.. np:function:: DMRAPFromMoments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`H = DMRAPFromMoments(moms, Nm)`
        * - Mathematica:
          - :code:`H = DMRAPFromMoments[moms, Nm]`
        * - Python/Numpy:
          - :code:`H = DMRAPFromMoments(moms, Nm)`

    Creates a discrete marked rational arrival process that
    has the same marginal and lag-1 joint moments as given 
    (see [1]_).

    Parameters
    ----------
    moms : vector of doubles
        The list of marginal moments. To obtain a discrete 
        marked rational process of order M, 2*M-1 marginal 
        moments are required.
    Nm : list of matrices, shape (M,M)
        The list of lag-1 joint moment matrices. The 
        length of the list determines K, the number of arrival 
        types of the discrete rational process.
    
    Returns
    -------
    H : list of matrices, shape (M,M)
        The H0, H1, ..., HK matrices of the discrete marked
        rational process

    References
    ----------
    .. [1] Andras Horvath, Gabor Horvath, Miklos Telek, "A 
           traffic based decomposition of two-class queueing
           networks with priority service," Computer Networks 
           53:(8) pp. 1235-1248. (2009)

    Examples
    ========
    For Matlab:

    >>> G0 = [0.34, 0, 0; 0.06, 0.05, 0.03; 0.11, 0.13, 0];
    >>> G1 = [0.3, 0, 0; 0.16, 0.18, 0.05; 0.15, 0.04, 0.09];
    >>> G2 = [0, 0.01, 0; 0.1, 0.07, 0.08; 0.13, 0.12, 0.13];
    >>> G3 = [0.35, 0, 0; 0, 0.18, 0.04; 0.06, 0.03, 0.01];
    >>> G = {G0, G1, G2, G3};
    >>> moms = MarginalMomentsFromDMRAP(G, 5);
    >>> disp(moms);
           1.5037       3.0278       8.4243       31.097       143.88
    >>> Nm = LagkJointMomentsFromDMRAP(G, 2, 1);
    >>> [Nm1, Nm2, Nm3] = Nm{:};
    >>> disp(Nm1);
          0.45395      0.68525       1.3856
          0.68283       1.0318       2.0887
           1.3755       2.0807       4.2171
    >>> disp(Nm2);
         0.026281      0.03323     0.053055
         0.035051     0.043866     0.068917
         0.060653      0.07477      0.11464
    >>> disp(Nm3);
          0.51977      0.78522       1.5891
          0.78582       1.1881       2.4067
           1.5917       2.4087       4.8838
    >>> H = DMRAPFromMoments(moms, Nm);
    >>> disp(H{1});
         0.067795      0.67922     -0.42509
        -0.018716    0.0035507      0.35386
        -0.019902     0.039925      0.31865
    >>> disp(H{2});
          0.19512       1.6571      -1.5453
       -0.0081004      0.18194      0.12658
       -0.0076905      0.11526      0.19294
    >>> disp(H{3});
          0.25787       3.8341      -4.0555
          0.12827      -1.3508       1.2347
          0.12666      -1.4076       1.2929
    >>> disp(H{4});
          0.17441     -0.11574      0.27606
        -0.014006    -0.013671      0.37636
        -0.012887      -0.0175      0.37926
    >>> BuToolsCheckPrecision = 10.^-10;
    >>> rmoms = MarginalMomentsFromDMRAP(H, 5);
    >>> disp(rmoms);
           1.5037       3.0278       8.4243       31.097       143.88
    >>> rNm = LagkJointMomentsFromDMRAP(H, 2, 1);
    >>> [rNm1, rNm2, rNm3] = rNm{:};
    >>> disp(rNm1);
          0.45395      0.68525       1.3856
          0.68283       1.0318       2.0887
           1.3755       2.0807       4.2171
    >>> disp(rNm2);
         0.026281      0.03323     0.053055
         0.035051     0.043866     0.068917
         0.060653      0.07477      0.11464
    >>> disp(rNm3);
          0.51977      0.78522       1.5891
          0.78582       1.1881       2.4067
           1.5917       2.4087       4.8838
    >>> disp(norm(moms-rmoms));
       3.0056e-11

    For Mathematica:

    
    For Python/Numpy:

    >>> G0 = ml.matrix([[0.34, 0, 0],[0.06, 0.05, 0.03],[0.11, 0.13, 0]])
    >>> G1 = ml.matrix([[0.3, 0, 0],[0.16, 0.18, 0.05],[0.15, 0.04, 0.09]])
    >>> G2 = ml.matrix([[0, 0.01, 0],[0.1, 0.07, 0.08],[0.13, 0.12, 0.13]])
    >>> G3 = ml.matrix([[0.35, 0, 0],[0, 0.18, 0.04],[0.06, 0.03, 0.01]])
    >>> G = [G0, G1, G2, G3]
    >>> moms = MarginalMomentsFromDMRAP(G, 5)
    >>> print(moms)
    [1.503697331491185, 3.0278418573508938, 8.424305390832199, 31.097386717744087, 143.88041840101141]
    >>> Nm = LagkJointMomentsFromDMRAP(G, 2, 1)
    >>> Nm1, Nm2, Nm3 = Nm
    >>> print(Nm1)
    [[ 0.45395  0.68525  1.38564]
     [ 0.68283  1.03179  2.08868]
     [ 1.37545  2.08069  4.21706]]
    >>> print(Nm2)
    [[ 0.02628  0.03323  0.05306]
     [ 0.03505  0.04387  0.06892]
     [ 0.06065  0.07477  0.11464]]
    >>> print(Nm3)
    [[ 0.51977  0.78522  1.58915]
     [ 0.78582  1.1881   2.40665]
     [ 1.59173  2.4087   4.88378]]
    >>> H = DMRAPFromMoments(moms, Nm)
    >>> print(H[0])
    [[ 0.0678   0.67922 -0.42509]
     [-0.01872  0.00355  0.35386]
     [-0.0199   0.03992  0.31865]]
    >>> print(H[1])
    [[ 0.19512  1.65707 -1.5453 ]
     [-0.0081   0.18194  0.12658]
     [-0.00769  0.11526  0.19294]]
    >>> print(H[2])
    [[ 0.25787  3.83412 -4.05553]
     [ 0.12827 -1.35077  1.2347 ]
     [ 0.12666 -1.40761  1.2929 ]]
    >>> print(H[3])
    [[ 0.17441 -0.11574  0.27606]
     [-0.01401 -0.01367  0.37636]
     [-0.01289 -0.0175   0.37926]]
    >>> butools.checkPrecision = 10.**-10
    >>> rmoms = MarginalMomentsFromDMRAP(H, 5)
    >>> print(rmoms)
    [1.5036973314922903, 3.0278418573555848, 8.4243053908503818, 31.09738671782155, 143.88041840139212]
    >>> rNm = LagkJointMomentsFromDMRAP(H, 2, 1)
    >>> rNm1, rNm2, rNm3 = rNm
    >>> print(rNm1)
    [[ 0.45395  0.68525  1.38564]
     [ 0.68283  1.03179  2.08868]
     [ 1.37545  2.08069  4.21706]]
    >>> print(rNm2)
    [[ 0.02628  0.03323  0.05306]
     [ 0.03505  0.04387  0.06892]
     [ 0.06065  0.07477  0.11464]]
    >>> print(rNm3)
    [[ 0.51977  0.78522  1.58915]
     [ 0.78582  1.1881   2.40665]
     [ 1.59173  2.4087   4.88378]]
    >>> print(la.norm(np.array(moms)-np.array(rmoms)))
    3.8896478315352875e-10

