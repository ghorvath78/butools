butools.fitting.RelativeEntropy
===============================

.. currentmodule:: butools.fitting

.. np:function:: RelativeEntropy

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`re = RelativeEntropy(p1, p2)`
        * - Mathematica:
          - :code:`re = RelativeEntropy[p1, p2]`
        * - Python/Numpy:
          - :code:`re = RelativeEntropy(p1, p2)`

    Returns the relative entropy (aka Kullback–Leibler
    divergence) of two vectors.

    Parameters
    ----------
    p1 : vector, length M
        The first vector
    p2 : vector, length M
        The second vector

    Returns
    -------
    re : double
        The relative entropy calculated as
        :math:`re=\sum_{i=1}^M p1_i |\log(p1_i/p2_i)|`

    Examples
    --------    
    For Matlab:
    
    >>> tr = dlmread('trace.txt');
    >>> trAcf = LagCorrelationsFromTrace(tr(1:10000), 10);
    >>> [D0,D1]=MAPFromTrace(tr(1:10000),5);
    >>> mAcf = LagCorrelationsFromMAP(D0, D1, 10, 1e-13);
    >>> RelativeEntropy (mAcf, trAcf)
          0.32132
    
    For Python/Numpy:
    
    >>> tr = np.loadtxt('trace.txt')
    >>> trAcf = LagCorrelationsFromTrace(tr[0:10000], 10)
    >>> D0, D1 = MAPFromTrace(tr[0:10000],5)
    >>> mAcf = LagCorrelationsFromMAP(D0, D1, 10, 1e-13)
    >>> print(RelativeEntropy (mAcf, trAcf))
    0.321362238531

