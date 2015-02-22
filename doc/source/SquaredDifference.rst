butools.fitting.SquaredDifference
=================================

.. currentmodule:: butools.fitting

.. np:function:: SquaredDifference

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`sd = SquaredDifference(p1, p2)`
        * - Mathematica:
          - :code:`sd = SquaredDifference[p1, p2]`
        * - Python/Numpy:
          - :code:`sd = SquaredDifference(p1, p2)`

    Returns the squared difference between two vectors.

    Parameters
    ----------
    p1 : vector, length M
        The first vector
    p2 : vector, length M
        The second vector

    Returns
    -------
    sd : double
        The squared difference calculated as
        :math:`sq=\sum_{i=1}^M (p1_i-p2_i)^2`

    Examples
    --------    
    For Matlab:
    
    >>> tr = dlmread('trace.txt');
    >>> trAcf = LagCorrelationsFromTrace(tr(1:10000), 10);
    >>> [D0,D1]=MAPFromTrace(tr(1:10000),5);
    >>> mAcf = LagCorrelationsFromMAP(D0, D1, 10, 1e-13);
    >>> SquaredDifference (mAcf, trAcf)
         0.023828
    
    For Python/Numpy:
    
    >>> tr = np.loadtxt('trace.txt')
    >>> trAcf = LagCorrelationsFromTrace(tr[0:10000], 10)
    >>> D0, D1 = MAPFromTrace(tr[0:10000],5)
    >>> mAcf = LagCorrelationsFromMAP(D0, D1, 10, 1e-13)
    >>> print(SquaredDifference (mAcf, trAcf))
    0.0238340444119

         
