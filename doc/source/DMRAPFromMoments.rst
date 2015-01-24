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

    For Python/Numpy:
    
    >>> moms = [[1.503697331491185, 3.0278418573508938, 8.424305390832199, 31.097386717744087, 143.88041840101141]
    >>> jmoms1 = ml.matrix([[ 0.45394628,  0.68525107,  1.38563721],[ 0.68282502,  1.03179221, 2.08867835], [ 1.37545474,  2.08069064,  4.21706252]])
    >>> jmoms2 = ml.matrix([[ 0.02628103,  0.03323045,  0.0530552 ],[ 0.03505148,  0.04386615,  0.06891662],[ 0.06065279,  0.07477024,  0.11464022]])
    >>> jmoms3 = ml.matrix([[ 0.51977269,  0.78521582,  1.58914945],[ 0.78582084,  1.18810084,  2.4066502 ],[ 1.59173433,  2.40869961,  4.88377952]])
    >>> H=DMRAPFromMoments(moms,Nm)
    >>> print(H[0])
    [[ 0.06779514  0.67922331 -0.42509134]
     [-0.01871617  0.00355068  0.35386442]
     [-0.01990152  0.0399245   0.31865418]]
    >>> print(H[1])
    [[ 0.19512078  1.65706952 -1.54529719]
     [-0.00810037  0.18194379  0.12657574]
     [-0.00769055  0.11525769  0.19293544]]
    >>> print(H[2])
    [[ 0.25786854  3.83411726 -4.05553027]
     [ 0.12827077 -1.35076999  1.23470023]
     [ 0.1266564  -1.40761393  1.29290145]]
    >>> print(H[3])
    [[ 0.17440749 -0.11574488  0.27606165]
     [-0.0140057  -0.01367117  0.37635777]
     [-0.01288707 -0.01750029  0.37926368]]
    >>> print(MarginalMomentsFromDMRAP(H,5,1e-10))
    [1.5036973314922903, 3.0278418573555848, 8.4243053908503818, 31.09738671782155, 143.88041840139212]
    >>> rNm1, rNm2, rNm3 = LagkJointMomentsFromDMRAP(H,2,1,1e-10)
    >>> print(rNm1)    
    [[ 0.45394628  0.68525107  1.38563721]
     [ 0.68282502  1.03179221  2.08867835]
     [ 1.37545474  2.08069064  4.21706252]]
    >>> print(rNm2)    
    [[ 0.02628103  0.03323045  0.0530552 ]
     [ 0.03505148  0.04386615  0.06891662]
     [ 0.06065279  0.07477024  0.11464022]]
    >>> print(rNm3)    
    [[ 0.51977269  0.78521582  1.58914945]
     [ 0.78582084  1.18810084  2.4066502 ]
     [ 1.59173433  2.40869961  4.88377952]]

