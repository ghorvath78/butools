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
    --------
    For Matlab:
    
    >>> moms = [1.5037, 3.02784, 8.42431, 31.0974, 143.88];
    >>> jmoms1 = [0.453946, 0.685251, 1.38564; 0.682825, 1.03179, 2.08868; 1.37545, 2.08069, 4.21706];
    >>> jmoms2 = [0.026281, 0.0332304, 0.0530552; 0.0350515, 0.0438661, 0.0689166; 0.0606528, 0.0747702, 0.11464];
    >>> jmoms3 = [0.519773, 0.785216, 1.58915; 0.785821, 1.1881, 2.40665; 1.59173, 2.4087, 4.88378];
    >>> H = DMRAPFromMoments(moms,{jmoms1,jmoms2,jmoms3});
    >>> H{1}
         0.067795      0.67922     -0.42509
        -0.018716    0.0035507      0.35386
        -0.019902     0.039925      0.31865
    >>> H{2}
          0.19512       1.6571      -1.5453
       -0.0081004      0.18194      0.12658
       -0.0076905      0.11526      0.19294
    >>> H{3}
          0.25787       3.8341      -4.0555
          0.12827      -1.3508       1.2347
          0.12666      -1.4076       1.2929
    >>> H{4}
          0.17441     -0.11574      0.27606
        -0.014006    -0.013671      0.37636
        -0.012887      -0.0175      0.37926
    >>> MarginalMomentsFromDMRAP(H,5,1e-11)
           1.5037       3.0278       8.4243       31.097       143.88
    >>> rjmoms = LagkJointMomentsFromDMRAP(H,2,1,1e-11);
    >>> rjmoms{1}
          0.45395      0.68525       1.3856
          0.68283       1.0318       2.0887
           1.3755       2.0807       4.2171
    >>> rjmoms{2}
         0.026281      0.03323     0.053055
         0.035051     0.043866     0.068917
         0.060653      0.07477      0.11464
    >>> rjmoms{3}
          0.51977      0.78522       1.5891
          0.78582       1.1881       2.4067
           1.5917       2.4087       4.8838

