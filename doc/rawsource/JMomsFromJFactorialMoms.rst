butools.moments.JMomsFromJFactorialMoms
=======================================

.. currentmodule:: butools.moments

.. np:function:: JMomsFromJFactorialMoms

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`jm = JMomsFromJFactorialMoms(jfm)`
        * - Mathematica:
          - :code:`jm = JMomsFromJFactorialMoms[jfm]`
        * - Python/Numpy:
          - :code:`jm = JMomsFromJFactorialMoms(jfm)`
    
    Returns the lag-1 joint raw moments given the lag-1 joint
    factorial moments.
    
    The lag-1 joint raw moments are: 
    :math:`m_{i,j}=E(\mathcal{X}^i\mathcal{Y}^j)`
    
    The factorial moments are: 
    :math:`f_{ij}=E(\mathcal{X}(\mathcal{X}-1)\cdots(\mathcal{X}-i+1)\mathcal{Y}(\mathcal{Y}-1)\cdots(\mathcal{Y}-j+1))`
       
    Parameters
    ----------
    jfm : matrix, shape (M,M)
        The matrix of joint factorial moments. The entry in 
        row i and column j is :math:`f_{i,j},i\geq 1,j\geq 1`.
        
    Returns
    -------
    jm : matrix, shape (M,M)
        The matrix of joint raw moments. The entry in row i
        and column j is :math:`m_{i,j},i\geq 1,j\geq 1`.
        
    References
    ----------
    http://en.wikipedia.org/wiki/Factorial_moment    

    Examples
    --------
    For Matlab:
    
    >>> jfm=JFactorialMomsFromJMoms([0.7 2 3 4; 5 6 7 8; 9 10 11 12])
    [0.7 1.3 -1.6 3.8; 4.3 -0.3 0.6 -1.8; -4.6 0.6 -1.2 3.6]
    >>> jm=JMomsFromJFactorialMoms(jfm)
    [0.7 2 3 4; 5 6 7 8; 9 10 11 12]
    
    For Mathematica:
    
    >>> jfm=JFactorialMomsFromJMoms[{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}}]
    {{1, 1, -1, 2}, {4, 0, 0, 0}, {-4, 0, 0, 0}}
    >>> jm=JMomsFromJFactorialMoms[jfm]
    {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}}
    
    For Python/Numpy:

    >>> jfm=JFactorialMomsFromJMoms([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
    >>> print(jfm)
    [[ 1.  1. -1.  2.]
     [ 4.  0.  0.  0.]
     [-4.  0.  0.  0.]]
    >>> jm=JMomsFromJFactorialMoms(jfm)
    >>> print(jm)
    [[  1.   2.   3.   4.]
     [  5.   6.   7.   8.]
     [  9.  10.  11.  12.]]
         
