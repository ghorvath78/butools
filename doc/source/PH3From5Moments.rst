butools.ph.PH3From5Moments
==========================

.. currentmodule:: butools.ph

.. np:function:: PH3From5Moments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = PH3From5Moments(moms)`
        * - Mathematica:
          - :code:`{alpha, A} = PH3From5Moments[moms]`
        * - Python/Numpy:
          - :code:`alpha, A = PH3From5Moments(moms)`

    Returns a PH(3) which has the same 5 moments as given.
    
    Parameters
    ----------
    moms : vector of doubles, length(5)
      The moments to match

    Returns
    -------
    alpha : vector, shape (1,3)
        The initial probability vector of the PH(3)
    A : matrix, shape (3,3)
        Transient generator matrix of the PH(3)

    Notes
    -----
    Raises an error if the moments are not feasible with
    a PH(3). Also note that the numerical behavior of the 
    procedure can be poor if the moments are close to the 
    boundary of the feasible region.

    References
    ----------
    .. [1] G. Horvath and M. Telek, "On the canonical 
           representation of phase type distributions," 
           Performance Evaluation, vol. 66, no. 8, pp. 
           396 - 409, 2009.

    Examples
    --------
    For Matlab:
    
    >>> moms = [0.20939, 0.10449, 0.089092, 0.11027, 0.17953];
    >>> [a,A]=PH3From5Moments(moms);
    >>> a
      0.57641       0.3324     0.091191   
    >>> A
      -10.102            0            0
       5.3797      -5.3797            0
            0       2.8798      -2.8798   
    >>> MomentsFromPH(a,A,5)
      0.20939      0.10449     0.089092      0.11027      0.17953

    >>> moms = [0.44865, 0.5496, 1.3298, 4.9428, 24.182];
    >>> [a,A]=PH3From5Moments(moms);
    >>> a
      0.96612     0.009501     0.024382
    >>> A
      -2.9245            0      0.12769
       2.8096      -2.8096            0
            0       1.1153      -1.1153
    >>> MomentsFromPH(a,A,5)
      0.44865       0.5496       1.3298       4.9428       24.182
      
    For Python/Numpy:
    
    >>> moms = [0.20939, 0.10449, 0.089092, 0.11027, 0.17953]
    >>> a,A=PH3From5Moments(moms)
    >>> print(a)
    [[ 0.57641145  0.33239799  0.09119055]]
    >>> print(A)
    [[-10.10242487   0.           0.        ]
     [  5.3797402   -5.3797402    0.        ]
     [  0.           2.87975038  -2.87975038]]
    >>> print(MomentsFromPH(a,A,5))
    [0.20938999999999752, 0.10448999999999731, 0.089091999999996285, 0.11026999999999371, 0.17952999999998739]

    >>> moms = [0.44865, 0.5496, 1.3298, 4.9428, 24.182]
    >>> a,A=PH3From5Moments(moms)
    >>> print(a)
    [[ 0.96611708  0.00950097  0.02438196]]
    >>> print(A)
    [[-2.92453876  0.          0.12768624]
     [ 2.80961128 -2.80961128  0.        ]
     [ 0.          1.11527694 -1.11527694]]
    >>> print(MomentsFromPH(a,A,5))
    [0.44865000000000199, 0.5496000000000052, 1.3298000000000183, 4.9428000000000774, 24.182000000000411]
              
