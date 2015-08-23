butools.dmap.LagkJointMomentsFromDRAP
=====================================

.. currentmodule:: butools.dmap

.. np:function:: LagkJointMomentsFromDRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`Nm = LagkJointMomentsFromDRAP(H0, H1, K, L, prec)`
        * - Mathematica:
          - :code:`Nm = LagkJointMomentsFromDRAP[H0, H1, K, L, prec]`
        * - Python/Numpy:
          - :code:`Nm = LagkJointMomentsFromDRAP(H0, H1, K, L, prec)`

    Returns the lag-L joint moments of a discrete rational 
    arrival process.

    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the discrete rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the discrete rational arrival process
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
    ========
    For Matlab:

    >>> H0 = [0, 0, 0.13; 0, 0.6, 0.18; 0.31, 0.26, 0.02];
    >>> H1 = [0, 1, -0.13; 0, 0.18, 0.04; 0.03, 0.09, 0.29];
    >>> Nm = LagkJointMomentsFromDRAP(H0, H1, 4, 1);
    >>> disp(length(Nm));
         5
    >>> moms = MarginalMomentsFromDRAP(H0, H1, 4);
    >>> disp(moms);
            3.207       16.898       130.77       1347.1
    >>> cjm = zeros(1,3);
    >>> for i=1:1:3
    >>>     Nx = LagkJointMomentsFromDRAP(H0, H1, 1, i);
    >>>     cjm(i) = (Nx(2, 2)-moms(1)^2)/(moms(2)-moms(1)^2);
    >>> end
    >>> disp(cjm);
         0.014303    0.0012424   7.5868e-06
    >>> corr = LagCorrelationsFromDRAP(H0, H1, 3);
    >>> disp(corr);
         0.014303    0.0012424   7.5868e-06

    For Mathematica:

    
    For Python/Numpy:

    >>> H0 = ml.matrix([[0, 0, 0.13],[0, 0.6, 0.18],[0.31, 0.26, 0.02]])
    >>> H1 = ml.matrix([[0, 1, -0.13],[0, 0.18, 0.04],[0.03, 0.09, 0.29]])
    >>> Nm = LagkJointMomentsFromDRAP(H0, H1, 4, 1)
    >>> print(Length(Nm))
    5
    >>> moms = MarginalMomentsFromDRAP(H0, H1, 4)
    >>> print(moms)
    [3.2070236684078202, 16.897636691953394, 130.77054574356021, 1347.0743893905096]
    >>> cjm = np.zeros(3)
    >>> for i in range(1,4,1):
    >>>     Nx = LagkJointMomentsFromDRAP(H0, H1, 1, i)
    >>>     cjm[i-1] = (Nx[1, 1]-moms[0]**2)/(moms[1]-moms[0]**2)
    >>> print(cjm)
    [  1.43030e-02   1.24240e-03   7.58676e-06]
    >>> corr = LagCorrelationsFromDRAP(H0, H1, 3)
    >>> print(corr)
    [  1.43030e-02   1.24240e-03   7.58676e-06]

