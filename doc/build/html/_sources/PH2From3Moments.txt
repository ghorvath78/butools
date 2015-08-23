butools.ph.PH2From3Moments
==========================

.. currentmodule:: butools.ph

.. np:function:: PH2From3Moments

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = PH2From3Moments(moms, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = PH2From3Moments[moms, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = PH2From3Moments(moms, prec)`

    Returns a PH(2) which has the same 3 moments as given.

    Parameters
    ----------
    moms : vector of doubles, length(3)
      The moments to match
    prec : double, optional
      Numerical precision, default value is 1e-14
    
    Returns
    -------
    alpha : matrix, shape (1,2)
      The initial probability vector of the PH(2)
    A : matrix, shape (2,2)
      Transient generator matrix of the PH(2)

    Notes
    -----
    Raises an error if the moments are not feasible with
    a PH(2).

    References
    ----------
    .. [1]  M. Telek and A. Heindl, "Moment bounds for acyclic 
            discrete and continuous phase-type distributions of
            second order," in In Proc. of UK Performance 
            Evaluation Workshop, UKPEW, 2002"

    Examples
    ========
    For Matlab:

    >>> moms = [10.0, 160.0, 3500.0];
    >>> [a, A] = PH2From3Moments(moms);
    >>> disp(a);
           0.8702       0.1298
    >>> disp(A);
         -0.15576      0.15576
                0     -0.22659
    >>> phmoms = MomentsFromPH(a, A, 3);
    >>> disp(phmoms);
              10         160        3500
    >>> moms = [10.0, 260.0, 13500.0];
    >>> [a, A] = PH2From3Moments(moms);
    >>> disp(a);
         0.090975      0.90902
    >>> disp(A);
        -0.041955     0.041955
                0     -0.12769
    >>> phmoms = MomentsFromPH(a, A, 3);
    >>> disp(phmoms);
               10          260        13500

    For Mathematica:

    
    For Python/Numpy:

    >>> moms = [10.0, 160.0, 3500.0]
    >>> a, A = PH2From3Moments(moms)
    >>> print(a)
    [[ 0.8702  0.1298]]
    >>> print(A)
    [[-0.15576  0.15576]
     [ 0.      -0.22659]]
    >>> phmoms = MomentsFromPH(a, A, 3)
    >>> print(phmoms)
    [10.0, 160.0, 3500.0]
    >>> moms = [10.0, 260.0, 13500.0]
    >>> a, A = PH2From3Moments(moms)
    >>> print(a)
    [[ 0.09098  0.90902]]
    >>> print(A)
    [[-0.04195  0.04195]
     [ 0.      -0.12769]]
    >>> phmoms = MomentsFromPH(a, A, 3)
    >>> print(phmoms)
    [10.0, 260.00000000000006, 13500.000000000005]

