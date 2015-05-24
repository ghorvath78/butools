butools.map.MarginalDistributionFromMAP
=======================================

.. currentmodule:: butools.map

.. np:function:: MarginalDistributionFromMAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = MarginalDistributionFromMAP(D0, D1, precision)`
        * - Mathematica:
          - :code:`{alpha, A} = MarginalDistributionFromMAP[D0, D1, precision]`
        * - Python/Numpy:
          - :code:`alpha, A = MarginalDistributionFromMAP(D0, D1, precision)`

    Returns the phase type distributed marginal distribution
    of a Markovian arrival process.

    Parameters
    ----------
    D0 : matrix, shape (M,M)
        The D0 matrix of the Markovian arrival process
    D1 : matrix, shape (M,M)
        The D1 matrix of the Markovian arrival process
    precision : double, optional
        Numerical precision for checking if the input is valid.
        The default value is 1e-14

    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial probability vector of the phase type 
        distributed marginal distribution
    A : matrix, shape (M,M)
        The transient generator of the phase type distributed
        marginal distribution    

    Examples
    --------
    For Matlab:
    
    >>> D0=[-0.17 0 0 0.07; 0.01 -0.78 0.03 0.08; 0.22 0.17 -1.1 0.02; 0.04 0.12 0 -0.42];
    >>> D1=[0 0.06 0 0.04; 0.04 0.19 0.21 0.22; 0.22 0.13 0.15 0.19; 0.05 0 0.17 0.04];
    >>> [a,A]=MarginalDistributionFromMAP(D0,D1);
    >>> a
          0.14438      0.23571      0.33794      0.28196
    >>> A
            -0.17            0            0         0.07
             0.01        -0.78         0.03         0.08
             0.22         0.17         -1.1         0.02
             0.04         0.12            0        -0.42
    >>> MomentsFromPH(a,A)
           3.4433        34.03       592.08        14548    4.559e+05    1.727e+07   7.6526e+08
    >>> x = (0:0.1:10);
    >>> plot(x,PdfFromPH(a,A,x));

    For Python/Numpy:
    
    >>> D0=ml.matrix([[-0.17, 0, 0, 0.07],[0.01, -0.78, 0.03, 0.08],[0.22, 0.17, -1.1, 0.02],[0.04, 0.12, 0, -0.42]])
    >>> D1=ml.matrix([[0, 0.06, 0, 0.04],[0.04, 0.19, 0.21, 0.22],[0.22, 0.13, 0.15, 0.19],[0.05, 0, 0.17, 0.04]])
    >>> [a,A]=MarginalDistributionFromMAP(D0,D1)
    >>> print(a)
    [[ 0.14438297  0.2357099   0.33794227  0.28196486]]    
    >>> print(A)
    [[-0.17  0.    0.    0.07]
     [ 0.01 -0.78  0.03  0.08]
     [ 0.22  0.17 -1.1   0.02]
     [ 0.04  0.12  0.   -0.42]]
    >>> print(MomentsFromPH(a,A))
    [3.4433473754205295, 34.030353434960716, 592.07698593172893, 14547.729458258538, 455900.27697389701, 17269527.005202383, 765259717.27887475]
    >>> x = np.arange(0,10,0.1)
    >>> plt.plot(x,PdfFromPH(a,A,x))

