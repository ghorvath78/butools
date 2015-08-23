butools.dph.CanonicalFromDPH3
=============================

.. currentmodule:: butools.dph

.. np:function:: CanonicalFromDPH3

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[beta, B] = CanonicalFromDPH3(alpha, A, prec)`
        * - Mathematica:
          - :code:`{beta, B} = CanonicalFromDPH3[alpha, A, prec]`
        * - Python/Numpy:
          - :code:`beta, B = CanonicalFromDPH3(alpha, A, prec)`

    Returns the canonical form of an order-3 discrete phase-type 
    distribution.

    Parameters
    ----------
    alpha : matrix, shape (1,3)
        Initial vector of the discrete phase-type distribution
    A : matrix, shape (3,3)
        Transition probability matrix of the discrete phase-type
        distribution
    prec : double, optional
      Numerical precision for checking the input, default value
      is 1e-14

    Returns
    -------
    beta : matrix, shape (1,3)
      The initial probability vector of the canonical form
    B : matrix, shape (3,3)
      Transition probability matrix of the canonical form

    Examples
    ========
    For Matlab:

    >>> a = [0.46,0.22,0.32];
    >>> A = [0.67, 0.01, 0.12; 0.06, 0.45, 0.15; 0.18, 0.43, 0.32];
    >>> [b, B] = CanonicalFromDPH3(a, A);
    >>> disp(b);
          0.21239      0.37004      0.41757
    >>> disp(B);
          0.10918            0            0
          0.45654      0.54346            0
                0      0.21265      0.78735
    >>> ev = eig(A);
    >>> disp(ev);
          0.78735
          0.54346
          0.10918
    >>> flag = CheckDPHRepresentation(b, B);
    >>> disp(flag);
         1
    >>> Cm = SimilarityMatrix(A, B);
    >>> err1 = norm(A*Cm-Cm*B);
    >>> err2 = norm(a*Cm-b);
    >>> disp(max(err1, err2));
       2.3226e-13
    >>> a = [0.76,0.12,0.12];
    >>> A = [0.31, 0., 0.; 0.98, 0., 0.02; 0.88, 0.04, 0.08];
    >>> [b, B] = CanonicalFromDPH3(a, A);
    >>> disp(b);
          0.13074       0.3444      0.52486
    >>> disp(B);
             0.31         0.69            0
                0         0.08         0.92
                0   0.00086957            0
    >>> ev = eig(A);
    >>> disp(ev);
          0.08899
       -0.0089898
             0.31
    >>> flag = CheckDPHRepresentation(b, B);
    >>> disp(flag);
         1
    >>> Cm = SimilarityMatrix(A, B);
    >>> err1 = norm(A*Cm-Cm*B);
    >>> err2 = norm(a*Cm-b);
    >>> disp(max(err1, err2));
        4.371e-16
    >>> a = [0.67,0.07,0.26];
    >>> A = [0.31, 0., 0.; 0.98, 0., 0.02; 0.88, 0.04, 0.08];
    >>> [b, B] = CanonicalFromDPH3(a, A);
    >>> disp(b);
          0.15814      0.37915       0.4627
    >>> disp(B);
             0.31         0.69            0
                0         0.08         0.92
                0   0.00086957            0
    >>> ev = eig(A);
    >>> disp(ev);
          0.08899
       -0.0089898
             0.31
    >>> flag = CheckDPHRepresentation(b, B);
    >>> disp(flag);
         1
    >>> Cm = SimilarityMatrix(A, B);
    >>> err1 = norm(A*Cm-Cm*B);
    >>> err2 = norm(a*Cm-b);
    >>> disp(max(err1, err2));
       4.2151e-16
    >>> a = [0.78,0.04,0.18];
    >>> A = [0.06, 0.25, 0.31; 0.45, 0.18, 0.33; 0.98, 0, 0.01];
    >>> [b, B] = CanonicalFromDPH3(a, A);
    >>> disp(b);
          0.43828      0.23849      0.32323
    >>> disp(B);
             0.25         0.75            0
          0.53747            0      0.46253
         0.072496            0            0
    >>> ev = eig(A);
    >>> disp(ev);
          0.79606
         -0.48028
        -0.065779
    >>> flag = CheckDPHRepresentation(b, B);
    >>> disp(flag);
         1
    >>> Cm = SimilarityMatrix(A, B);
    >>> err1 = norm(A*Cm-Cm*B);
    >>> err2 = norm(a*Cm-b);
    >>> disp(max(err1, err2));
       1.4697e-15

    For Mathematica:

    
    For Python/Numpy:

    >>> a = ml.matrix([[0.46,0.22,0.32]])
    >>> A = ml.matrix([[0.67, 0.01, 0.12],[0.06, 0.45, 0.15],[0.18, 0.43, 0.32]])
    >>> b, B = CanonicalFromDPH3(a, A)
    >>> print(b)
    [[ 0.21239  0.37004  0.41757]]
    >>> print(B)
    [[ 0.10918  0.       0.     ]
     [ 0.45654  0.54346  0.     ]
     [ 0.       0.21265  0.78735]]
    >>> ev = la.eigvals(A)
    >>> print(ev)
    [ 0.78735+0.j  0.54346+0.j  0.10918+0.j]
    >>> flag = CheckDPHRepresentation(b, B)
    >>> print(flag)
    True
    >>> Cm = SimilarityMatrix(A, B)
    >>> err1 = la.norm(A*Cm-Cm*B)
    >>> err2 = la.norm(a*Cm-b)
    >>> print(np.max(err1, err2))
    6.85597820868e-13
    >>> a = ml.matrix([[0.76,0.12,0.12]])
    >>> A = ml.matrix([[0.31, 0., 0.],[0.98, 0., 0.02],[0.88, 0.04, 0.08]])
    >>> b, B = CanonicalFromDPH3(a, A)
    >>> print(b)
    [[ 0.13074  0.3444   0.52486]]
    >>> print(B)
    [[  3.10000e-01   6.90000e-01   0.00000e+00]
     [  0.00000e+00   8.00000e-02   9.20000e-01]
     [  0.00000e+00   8.69565e-04   0.00000e+00]]
    >>> ev = la.eigvals(A)
    >>> print(ev)
    [ 0.08899+0.j -0.00899+0.j  0.31000+0.j]
    >>> flag = CheckDPHRepresentation(b, B)
    >>> print(flag)
    True
    >>> Cm = SimilarityMatrix(A, B)
    >>> err1 = la.norm(A*Cm-Cm*B)
    >>> err2 = la.norm(a*Cm-b)
    >>> print(np.max(err1, err2))
    3.70825293136e-16
    >>> a = ml.matrix([[0.67,0.07,0.26]])
    >>> A = ml.matrix([[0.31, 0., 0.],[0.98, 0., 0.02],[0.88, 0.04, 0.08]])
    >>> b, B = CanonicalFromDPH3(a, A)
    >>> print(b)
    [[ 0.15814  0.37915  0.4627 ]]
    >>> print(B)
    [[  3.10000e-01   6.90000e-01   0.00000e+00]
     [  0.00000e+00   8.00000e-02   9.20000e-01]
     [  0.00000e+00   8.69565e-04   0.00000e+00]]
    >>> ev = la.eigvals(A)
    >>> print(ev)
    [ 0.08899+0.j -0.00899+0.j  0.31000+0.j]
    >>> flag = CheckDPHRepresentation(b, B)
    >>> print(flag)
    True
    >>> Cm = SimilarityMatrix(A, B)
    >>> err1 = la.norm(A*Cm-Cm*B)
    >>> err2 = la.norm(a*Cm-b)
    >>> print(np.max(err1, err2))
    3.70825293136e-16
    >>> a = ml.matrix([[0.78,0.04,0.18]])
    >>> A = ml.matrix([[0.06, 0.25, 0.31],[0.45, 0.18, 0.33],[0.98, 0, 0.01]])
    >>> b, B = CanonicalFromDPH3(a, A)
    >>> print(b)
    [[ 0.43828  0.23849  0.32323]]
    >>> print(B)
    [[ 0.25     0.75     0.     ]
     [ 0.53747  0.       0.46253]
     [ 0.0725   0.       0.     ]]
    >>> ev = la.eigvals(A)
    >>> print(ev)
    [ 0.79606+0.j -0.48028+0.j -0.06578+0.j]
    >>> flag = CheckDPHRepresentation(b, B)
    >>> print(flag)
    True
    >>> Cm = SimilarityMatrix(A, B)
    >>> err1 = la.norm(A*Cm-Cm*B)
    >>> err2 = la.norm(a*Cm-b)
    >>> print(np.max(err1, err2))
    8.57844876274e-16

