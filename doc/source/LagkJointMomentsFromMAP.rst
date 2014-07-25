butools.map.LagkJointMomentsFromMAP
===================================

.. currentmodule:: butools.map

.. np:function:: LagkJointMomentsFromMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Nm = LagkJointMomentsFromMAP(D0, D1, K, L, prec)`
        * - Mathematica:
          - :code:`Nm = LagkJointMomentsFromMAP[D0, D1, K, L, prec]`
        * - Python/Numpy:
          - :code:`Nm = LagkJointMomentsFromMAP(D0, D1, K, L, prec)`

    Returns the lag-L joint moments of a Markovian arrival
    process.

    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process
    K : int, optional
        The dimension of the matrix of joint moments to 
        compute. If K=0, the MxM joint moments will be 
        computed. The default value is 0
    L : int, optional
        The lag at which the joint moments are computed.
        The default value is 1
    prec : double, optional
        Numerical precision to check if the input is valid. 
        The default value is 1e-14

    Returns
    -------
    Nm : matrix, shape(K+1,K+1)
        The matrix containing the lag-L joint moments, 
        starting from moment 0.

    Examples
    --------
    For Matlab:
    
    >>> D0=[-5 0 1 1; 1 -8 1 0; 1 0 -4 1; 1 2 3 -9];
    >>> D1=[0 1 0 2; 2 1 3 0; 0 0 1 1; 1 1 0 1];
    >>> LagkJointMomentsFromMAP(D0,D1,4,1)
                1      0.34247      0.25054      0.28271      0.42984
          0.34247       0.1173     0.085789     0.096807      0.14721
          0.25054       0.0857     0.062633      0.07066      0.10744
          0.28271     0.096627     0.070589     0.079623      0.12107
          0.42984      0.14686      0.10727      0.12099      0.18396
    >>> MarginalMomentsFromMAP(D0,D1,4)
          0.34247      0.25054      0.28271      0.42984    

