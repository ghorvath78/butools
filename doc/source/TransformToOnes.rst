butools.reptrans.TransformToOnes
================================

.. currentmodule:: butools.reptrans

.. np:function:: TransformToOnes

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`B = TransformToOnes(clovec)`
        * - Mathematica:
          - :code:`B = TransformToOnes[clovec]`
        * - Python/Numpy:
          - :code:`B = TransformToOnes(clovec)`
    
    Returns the similarity transformation matrix that converts 
    the given column vector to a vector of ones. It works even
    if it has zero entries.
    
    Parameters
    ----------
    clovec : column vector, shape(M,1)
        The original closing vector
        
    Returns
    -------
    B : matrix, shape(M,M)
        The matrix by which :math:`B\cdot clovec = \mathbf{1}` holds

    Examples
    --------
    For Matlab:
    
    >>> clovec = [0.0, 0.3, -1.5, 0.0]';
    >>> B = TransformToOnes (clovec)
            0       3.3333            0            0
       3.3333       3.3333            0            0
       3.3333       3.3333            0       3.3333
     -0.83333     -0.83333     -0.83333     -0.83333
    >>> B*clovec
       1
       1
       1
       1

