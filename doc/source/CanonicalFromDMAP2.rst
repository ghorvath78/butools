butools.dmap.CanonicalFromDMAP2
===============================

.. currentmodule:: butools.dmap

.. np:function:: CanonicalFromDMAP2

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[G0, G1] = CanonicalFromDMAP2(D0, D1, prec)`
        * - Mathematica:
          - :code:`{G0, G1} = CanonicalFromDMAP2[D0, D1, prec]`
        * - Python/Numpy:
          - :code:`G0, G1 = CanonicalFromDMAP2(D0, D1, prec)`

    Returns the canonical form of an order-2 discrete Markovian
    arrival process.

    Parameters
    ----------
    D0 : matrix, shape (2,2)
        The D0 matrix of the DMAP(2)
    D1 : matrix, shape (2,2)
        The D1 matrix of the DMAP(2)
    prec : double, optional
        Numerical precision to check the input, default 
        value is 1e-14

    Returns
    -------
    G0 : matrix, shape (1,2)
        The D0 matrix of the canonical DMAP(2)
    G1 : matrix, shape (2,2)
        The D1 matrix of the canonical DMAP(2)

    Examples
    --------
    For Matlab:
    
    >>> D0=[0.26 0.28; 0.35 0.23]
    >>> D1=[0.28 0.18; 0.14 0.28]
    >>> [H0,H1]=CanonicalFromDMAP2(D0,D1);
    >>> H0
             0.49      0.38875
         0.098265            0
    >>> H1
          0.12125            0
          0.46299      0.43875
    >>> C=SimilarityMatrix(H0,D0);
    >>> dissimilarity = norm(H0*C-C*D0) + norm(H1*C-C*D1)
        1.551e-14

