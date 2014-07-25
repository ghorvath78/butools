butools.dmap.DMAPFromDRAP
=========================

.. currentmodule:: butools.dmap

.. np:function:: DMAPFromDRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[D0, D1] = DMAPFromDRAP(H0, H1, precision)`
        * - Mathematica:
          - :code:`{D0, D1} = DMAPFromDRAP[H0, H1, precision]`
        * - Python/Numpy:
          - :code:`D0, D1 = DMAPFromDRAP(H0, H1, precision)`

    Obtains a Markovian representation of a discrete rational
    arrival process of the same size, if possible, using the
    procedure published in [1]_.

    Parameters
    ----------
    H0 : matrix, shape (M,M)
        The H0 matrix of the discrete rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the discrete rational arrival process
    precision : double, optional
        A representation is considered to be a Markovian one
        if it is closer to it than this precision

    Returns
    -------
    D0 : matrix, shape (M,M)
        The D0 matrix of the discrete Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the discrete Markovian arrival process

    References
    ----------
    .. [1] G Horvath, M Telek, "A minimal representation of 
           Markov arrival processes and a moments matching 
           method," Performance Evaluation 64:(9-12) pp. 
           1153-1168. (2007)       

    Examples
    --------
    For Matlab:
    
    >>> H0=[0 0 0.13; 0 0.6 0.18; 0.31 0.26 0.02];
    >>> H1=[0 1 -0.13; 0 0.18 0.04; 0.03 0.09 0.29];
    >>> CheckDMAPRepresentation(H0,H1)
    CheckDMAPRepresentation: One of the matrices has negative element!
         0
    >>> [D0,D1]=DMAPFromDRAP(H0,H1);
    >>> D0
         0.051945     0.066321      0.12704
         0.011717      0.56745      0.29444
          0.41438      0.17501   0.00060547
    >>> D1
         0.085648      0.64664       0.0224
        0.0054434     0.089137     0.031816
          0.04656     0.068225      0.29521
    >>> CheckDMAPRepresentation(D0,D1)
         1
    >>> error = norm(LagkJointMomentsFromDRAP(H0,H1,4,1)-LagkJointMomentsFromDMAP(D0,D1,4,1))
       1.1971e-08
    
