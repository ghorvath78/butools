butools.dmap.DMMAPFromDMRAP
===========================

.. currentmodule:: butools.dmap

.. np:function:: DMMAPFromDMRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`D = DMMAPFromDMRAP(H, precision)`
        * - Mathematica:
          - :code:`D = DMMAPFromDMRAP[H, precision]`
        * - Python/Numpy:
          - :code:`D = DMMAPFromDMRAP(H, precision)`

    Obtains a Markovian representation of a discrete rational
    arrival process of the same size, if possible, using the
    procedure published in [1]_.

    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the DMRAP to transform
    precision : double, optional
        A representation is considered to be a Markovian one
        if it is closer to it than this precision

    Returns
    -------
    D : list/cell of matrices of shape(M,M), length(N)
        The D0...DN matrices of the DMMAP (if found)

    References
    ----------
    .. [1] András Horváth, Gábor Horváth, Miklós Telek, "A 
           traffic based decomposition of two-class queueing 
           networks with priority service". COMPUTER NETWORKS 
           53:(8) pp. 1235-1248. (2009)

    Examples
    --------
    For Matlab:

    >>> H0=[0.15 0.2 0.18; -0.20 0.17 0.22; 0.19 0.15 0.16];
    >>> H1=[0.01 0.08 0.16; 0.02 0.2 0.07; 0.02 0.15 0.17];
    >>> H2=[0.14 0.07 0.01; 0.19 0.02 0.31; 0.06 0.1 0];
    >>> H = {H0,H1,H2};
    >>> CheckDMMAPRepresentation(H)
    CheckProbMatrix: the matrix has negative element (precision: 1e-14)!
    CheckDMAPRepresentation: D0 isn't a transient probability matrix!
         0
    >>> D=DMMAPFromDMRAP(H);
    >>> D{1}
          0.12467      0.28836      0.11942
        6.707e-05       0.1772     0.069012
          0.14879      0.12427      0.17813
    >>> D{2}
         0.029906    6.864e-05      0.12993
         0.065653      0.26951   0.00057292
         0.074333      0.21952     0.080583
    >>> D{3}
          0.15173    0.0059825      0.14994
           0.1319    0.0079169      0.27817
          0.04554      0.12849   0.00035384
    >>> CheckDMMAPRepresentation(D)
         1
    >>> jmom=LagkJointMomentsFromDMRAP(H,3,1);
    >>> rjmom=LagkJointMomentsFromDMMAP(D,3,1);
    >>> error = norm(rjmom{1}-jmom{1}) + norm(rjmom{2}-jmom{2})
        1.748e-13
            
