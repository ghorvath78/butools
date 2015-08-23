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
    ========
    For Matlab:

    >>> a = [0.1,0.9,0];
    >>> A = [-6.2, 2, 0; 2, -9, 1; 1, 0, -3];
    >>> moms = MomentsFromPH(a, A);
    >>> disp(moms);
          0.20939      0.10449     0.089092      0.11027      0.17953
    >>> [a, A] = PH3From5Moments(moms);
    >>> disp(a);
          0.58305      0.32736     0.089589
    >>> disp(A);
          -9.9819            0            0
           5.3405      -5.3405            0
                0       2.8776      -2.8776
    >>> phmoms = MomentsFromME(a, A, 5);
    >>> disp(phmoms);
          0.20939      0.10449     0.089092      0.11027      0.17953
    >>> a = [0.2,0.3,0.5];
    >>> A = [-1, 0, 0; 0, -3, 0.5; 0, -0.5, -3];
    >>> moms = MomentsFromME(a, A);
    >>> disp(moms);
          0.44865       0.5496       1.3298       4.9428       24.182
    >>> [a, A] = PH3From5Moments(moms);
    >>> disp(a);
          0.94865     0.036778     0.014574
    >>> disp(A);
               -3            0      0.15385
            2.866       -2.866            0
                0        1.134       -1.134
    >>> phmoms = MomentsFromME(a, A, 5);
    >>> disp(phmoms);
          0.44865       0.5496       1.3298       4.9428       24.182

    For Mathematica:

    
    For Python/Numpy:

    >>> a = ml.matrix([[0.1,0.9,0]])
    >>> A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
    >>> moms = MomentsFromPH(a, A)
    >>> print(moms)
    [0.20938722294654497, 0.10448912014333091, 0.089091500391907288, 0.11026674096545433, 0.17953027324720897]
    >>> a, A = PH3From5Moments(moms)
    >>> print(a)
    [[ 0.58305  0.32736  0.08959]]
    >>> print(A)
    [[-9.98192  0.       0.     ]
     [ 5.34047 -5.34047  0.     ]
     [ 0.       2.87761 -2.87761]]
    >>> phmoms = MomentsFromME(a, A, 5)
    >>> print(phmoms)
    [0.20938722294654719, 0.10448912014333331, 0.089091500391910647, 0.11026674096546003, 0.17953027324722048]
    >>> a = ml.matrix([[0.2,0.3,0.5]])
    >>> A = ml.matrix([[-1, 0, 0],[0, -3, 0.5],[0, -0.5, -3]])
    >>> moms = MomentsFromME(a, A)
    >>> print(moms)
    [0.44864864864864862, 0.54959824689554426, 1.3298244921327464, 4.9427685241556087, 24.182331446704147]
    >>> a, A = PH3From5Moments(moms)
    >>> print(a)
    [[ 0.94865  0.03678  0.01457]]
    >>> print(A)
    [[-3.       0.       0.15385]
     [ 2.86603 -2.86603  0.     ]
     [ 0.       1.13397 -1.13397]]
    >>> phmoms = MomentsFromME(a, A, 5)
    >>> print(phmoms)
    [0.44864864864864928, 0.54959824689554582, 1.3298244921327509, 4.9427685241556283, 24.18233144670425]

