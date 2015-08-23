butools.map.MarginalMomentsFromMRAP
===================================

.. currentmodule:: butools.map

.. np:function:: MarginalMomentsFromMRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`moms = MarginalMomentsFromMRAP(H, K, precision)`
        * - Mathematica:
          - :code:`moms = MarginalMomentsFromMRAP[H, K, precision]`
        * - Python/Numpy:
          - :code:`moms = MarginalMomentsFromMRAP(H, K, precision)`

    Returns the moments of the marginal distribution of a 
    marked rational arrival process.

    Parameters
    ----------
    H : list/cell of matrices of shape(M,M), length(N)
        The H0...HN matrices of the MRAP
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

    >>> x = 0.18;
    >>> H0 = [-5., 0.1+x, 0.9, 1.; 1., -8., 0.9, 0.1; 0.9, 0.1, -4., 1.; 1., 2., 3., -9.];
    >>> H1 = [0.1-x, 0.7, 0.1, 0.1; 0.1, 1., 1.8, 0.1; 0.1, 0.1, 0.1, 0.7; 0.7, 0.1, 0.1, 0.1];
    >>> H2 = [0.1, 0.1, 0.1, 1.7; 1.8, 0.1, 1., 0.1; 0.1, 0.1, 0.7, 0.1; 0.1, 1., 0.1, 0.8];
    >>> moms = MarginalMomentsFromMRAP({H0, H1, H2});
    >>> disp(moms);
      Columns 1 through 6
          0.33951      0.24583      0.27424      0.41206      0.77677       1.7594
      Column 7
           4.6515

    For Mathematica:

    
    For Python/Numpy:

    >>> x = 0.18
    >>> H0 = ml.matrix([[-5., 0.1+x, 0.9, 1.],[1., -8., 0.9, 0.1],[0.9, 0.1, -4., 1.],[1., 2., 3., -9.]])
    >>> H1 = ml.matrix([[0.1-x, 0.7, 0.1, 0.1],[0.1, 1., 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
    >>> H2 = ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1., 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1., 0.1, 0.8]])
    >>> moms = MarginalMomentsFromMRAP([H0, H1, H2])
    >>> print(moms)
    [0.33950747762450084, 0.24582557198236554, 0.27423742766051129, 0.41206018133500932, 0.7767718404933559, 1.7594286078546524, 4.6515347631617807]

