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
    ========
    For Matlab:

    >>> MM = [0.7, 2., 3., 4.; 5., 6., 7., 8.; 9., 10., 11., 12.];
    >>> JFmoms = JFactorialMomsFromJMoms(MM);
    >>> disp(JFmoms);
              0.7          1.3         -1.6          3.8
              4.3         -0.3          0.6         -1.8
             -4.6          0.6         -1.2          3.6
    >>> Jmoms = JMomsFromJFactorialMoms(JFmoms);
    >>> disp(Jmoms);
              0.7            2            3            4
                5            6            7            8
                9           10           11           12
    >>> err = norm(Jmoms-MM);
    >>> disp(err);
         0

    For Mathematica:

    >>> MM = {{0.7, 2., 3., 4.},{5., 6., 7., 8.},{9., 10., 11., 12.}};
    >>> JFmoms = JFactorialMomsFromJMoms[MM];
    >>> Print[JFmoms];
    {{0.7, 1.3, -1.5999999999999996, 3.8000000000000007},
     {4.3, -0.30000000000000004, 0.5999999999999996, -1.8000000000000007},
     {-4.6, 0.6000000000000001, -1.1999999999999993, 3.6000000000000014}}
    >>> Jmoms = JMomsFromJFactorialMoms[JFmoms];
    >>> Print[Jmoms];
    {{0.7, 2., 3., 4.},
     {5., 6., 7., 8.},
     {9., 10., 11., 12.}}
    >>> err = Norm[Jmoms-MM];
    >>> Print[err];
    0.

    For Python/Numpy:

    >>> MM = ml.matrix([[0.7, 2., 3., 4.],[5., 6., 7., 8.],[9., 10., 11., 12.]])
    >>> JFmoms = JFactorialMomsFromJMoms(MM)
    >>> print(JFmoms)
    [[ 0.7  1.3 -1.6  3.8]
     [ 4.3 -0.3  0.6 -1.8]
     [-4.6  0.6 -1.2  3.6]]
    >>> Jmoms = JMomsFromJFactorialMoms(JFmoms)
    >>> print(Jmoms)
    [[  0.7   2.    3.    4. ]
     [  5.    6.    7.    8. ]
     [  9.   10.   11.   12. ]]
    >>> err = la.norm(np.array(Jmoms)-MM)
    >>> print(err)
    0.0

