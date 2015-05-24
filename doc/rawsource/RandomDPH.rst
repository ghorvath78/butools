butools.dph.RandomDPH
=====================

.. currentmodule:: butools.dph

.. np:function:: RandomDPH

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = RandomDPH(order, mean, zeroEntries, maxTrials, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = RandomDPH[order, mean, zeroEntries, maxTrials, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = RandomDPH(order, mean, zeroEntries, maxTrials, prec)`

    Returns a random discrete phase-type distribution with a 
    given mean value.

    Parameters
    ----------
    order : int
        The size of the discrete phase-type distribution
    mean : double, optional
        The mean of the discrete phase-type distribution 
    zeroEntries : int, optional
        The number of zero entries in the initial vector, 
        generator matrix and closing vector
    maxTrials : int, optional
        The maximum number of trials to find a proper DPH 
        (that has an irreducible phase process and none of 
        its parameters is all-zero). The default value is 
        1000.
    prec : double, optional
        Numerical precision for checking the irreducibility.
        The default value is 1e-14.

    Returns
    -------
    alpha : vector, shape (1,M)
        The initial probability vector of the phase-type 
        distribution.
    A : matrix, shape (M,M)
        The transient generator matrix of the phase-type 
        distribution.

    Notes
    -----
    If the procedure fails, try to increase the 'maxTrials'
    parameter, or increase the mean value.

    Examples
    --------
    For Matlab:

    >>> [alpha,A]=RandomDPH(3,10,5);
    >>> alpha
           0.87199            0      0.12801
    >>> A
          0.73786      0.26214            0
                0      0.94274     0.057263
         0.053768      0.19917       0.7011
    >>> 1-sum(A,2)
                0
                0
         0.045964

    For Python/Numpy:
    
    >>> a,A=RandomDPH(3,10,5)
    >>> print(a)
    [[ 0.23191789  0.76808211  0.        ]]
    >>> print(A)
    [[ 0.37841386  0.23180007  0.13802   ]
     [ 0.36207828  0.63792172  0.        ]
     [ 0.40193095  0.          0.59806905]]
    >>> print(1-np.sum(A,1))
    [[ 0.25176607]
     [ 0.        ]
     [ 0.        ]]
 

