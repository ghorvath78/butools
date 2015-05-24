butools.map.RAPFromMomentsAndCorrelations
=========================================

.. currentmodule:: butools.map

.. np:function:: RAPFromMomentsAndCorrelations

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[H0, H1] = RAPFromMomentsAndCorrelations(moms, corr)`
        * - Mathematica:
          - :code:`{H0, H1} = RAPFromMomentsAndCorrelations[moms, corr]`
        * - Python/Numpy:
          - :code:`H0, H1 = RAPFromMomentsAndCorrelations(moms, corr)`

    Returns a rational arrival process that has the same moments
    and lag autocorrelation coefficients as given.

    Parameters
    ----------
    moms : vector of doubles
        The vector of marginal moments. To obtain a RAP of 
        size M, 2*M-1 moments are required.
    corr : vector of doubles
        The vector of lag autocorrelation coefficients. To 
        obtain a RAP of size M, 2*M-3 coefficients are needed.
    
    Returns
    -------
    H0 : matrix, shape (M,M)
        The H0 matrix of the rational arrival process
    H1 : matrix, shape (M,M)
        The H1 matrix of the rational arrival process

    Notes
    -----
    There is no guarantee that the returned matrices define
    a valid stochastic process. The joint densities may be
    negative.

    References
    ----------
    .. [1] Mitchell, Kenneth, and Appie van de Liefvoort. 
           "Approximation models of feed-forward G/G/1/N 
           queueing networks with correlated arrivals." 
           Performance Evaluation 51.2 (2003): 137-152.

    Examples
    --------
    For Matlab:
    
    >>> moms = [0.29774, 0.19284, 0.19448, 0.26597, 0.45833];
    >>> corr = [0.012394, 0.0027412, 0.00072384];
    >>> [G0,G1]=RAPFromMomentsAndCorrelations(mom,corr);
    >>> G0
          -8.9629       22.253      -18.544
         -0.99178       -4.667        2.331
          -1.2473       2.4279      -4.5701
    >>> G1
            2.203      -1.3184       4.3699
           1.2179       1.8219      0.28791
           1.0212      0.41715       1.9512
    >>> MarginalMomentsFromRAP(G0,G1)
          0.29774      0.19284      0.19448      0.26597      0.45833
    >>> LagCorrelationsFromRAP(G0,G1,3)
         0.012394
        0.0027412
       0.00072384
    
    For Python/Numpy:
    
    >>> moms = [0.29774127310061604, 0.19283643304803644, 0.19448147792730755, 0.26597325539245531, 0.45833053059627116]
    >>> corr = [0.012393574884970124, 0.0027411959690404088, 0.00072383642135696988]
    >>> [G0,G1]=RAPFromMomentsAndCorrelations(mom,corr)
    >>> print(G0)
    [[ -8.96289388  22.25257011 -18.54409809]
     [ -0.99178156  -4.66699225   2.33103341]
     [ -1.2472989    2.42791179  -4.57011387]]
    >>> print(G1)
    [[ 2.20274746 -1.3172514   4.36892581]
     [ 1.2179263   1.82172664  0.28808746]
     [ 1.02115416  0.41735382  1.950993  ]]
    >>> print(MarginalMomentsFromRAP(G0,G1))
    [0.29774127310061604, 0.19283643304803638, 0.19448147792730741, 0.26597325539245492, 0.45833053059627044]
    >>> print(LagCorrelationsFromRAP(G0,G1,3))
    [0.012393574884970265, 0.0027411959690409431, 0.00072383642135750309]

