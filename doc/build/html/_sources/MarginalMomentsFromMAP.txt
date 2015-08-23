butools.map.MarginalMomentsFromMAP
==================================

.. currentmodule:: butools.map

.. np:function:: MarginalMomentsFromMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`moms = MarginalMomentsFromMAP(D0, D1, K, precision)`
        * - Mathematica:
          - :code:`moms = MarginalMomentsFromMAP[D0, D1, K, precision]`
        * - Python/Numpy:
          - :code:`moms = MarginalMomentsFromMAP(D0, D1, K, precision)`

    Returns the moments of the marginal distribution of a 
    Markovian arrival process.

    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process
    K : int, optional
        Number of moments to compute. If K=0, 2*M-1 moments
        are computed. The default value is K=0.
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14

    Returns
    -------
    moms : row vector of doubles, length K
        The vector of moments.

    Examples
    ========
    For Matlab:

    >>> D0 = [-0.17, 0, 0, 0.07; 0.01, -0.78, 0.03, 0.08; 0.22, 0.17, -1.1, 0.02; 0.04, 0.12, 0, -0.42];
    >>> D1 = [0, 0.06, 0, 0.04; 0.04, 0.19, 0.21, 0.22; 0.22, 0.13, 0.15, 0.19; 0.05, 0, 0.17, 0.04];
    >>> moms = MarginalMomentsFromMAP(D0, D1);
    >>> disp(moms);
      Columns 1 through 6
           3.4433        34.03       592.08        14548    4.559e+05    1.727e+07
      Column 7
       7.6526e+08

    For Mathematica:

    
    For Python/Numpy:

    >>> D0 = ml.matrix([[-0.17, 0, 0, 0.07],[0.01, -0.78, 0.03, 0.08],[0.22, 0.17, -1.1, 0.02],[0.04, 0.12, 0, -0.42]])
    >>> D1 = ml.matrix([[0, 0.06, 0, 0.04],[0.04, 0.19, 0.21, 0.22],[0.22, 0.13, 0.15, 0.19],[0.05, 0, 0.17, 0.04]])
    >>> moms = MarginalMomentsFromMAP(D0, D1)
    >>> print(moms)
    [3.4433473754205295, 34.030353434960716, 592.07698593172893, 14547.729458258538, 455900.27697389701, 17269527.005202383, 765259717.27887475]

