butools.map.MinimalRepFromMRAP
==============================

.. currentmodule:: butools.map

.. np:function:: MinimalRepFromMRAP

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`D = MinimalRepFromMRAP(H, how, precision)`
        * - Mathematica:
          - :code:`H = MinimalRepFromMRAP[H, how, precision]`
        * - Python/Numpy:
          - :code:`H = MinimalRepFromMRAP(H, how, precision)`

    Returns the minimal representation of a marked rational
    arrival process.

    Parameters
    ----------
    H : list of matrices of shape (M,M)
        The list of H0, H1, ..., HK matrices of the marked
        rational arrival process
    how : {"obs", "cont", "obscont"}, optional        
        Determines how the representation is minimized. 
        "cont" means controllability, "obs" means 
        observability, "obscont" means that the rational arrival
        process is minimized in both respects. Default value 
        is "obscont".
    precision : double, optional
       Precision used by the Staircase algorithm. The default
       value is 1e-12.
    
    Returns
    -------
    D : list of matrices of shape (M,M)
        The D0, D1, ..., DK matrices of the minimal 
        representation

    References
    ----------
    .. [1] P. Buchholz, M. Telek, "On minimal representation of 
           rational arrival processes." Madrid Conference on 
           Qeueuing theory (MCQT), June 2010.

    Examples
    ========
    For Matlab:

    >>> D0 = [-5., 1., 0; 3., -3., 0; 1., 1., -5.];
    >>> D1 = [0, 0, 0.8; 0, 0, 0; 0.2, 0.2, 0.2];
    >>> D2 = [0, 0, 3.2; 0, 0, 0; 0.8, 0.8, 0.8];
    >>> Dm = {D0, D1, D2};
    >>> H = MinimalRepFromMRAP(Dm, 'cont');
    >>> disp(H{1});
        -5     1     0
         3    -3     0
         1     1    -5
    >>> disp(H{2});
                0            0          0.8
                0            0            0
              0.2          0.2          0.2
    >>> disp(H{3});
                0            0          3.2
                0            0            0
              0.8          0.8          0.8
    >>> Cm = SimilarityMatrix(H{1}, Dm{1});
    >>> err = norm(H{1}*Cm-Cm*Dm{1})+norm(H{2}*Cm-Cm*Dm{2})+norm(H{3}*Cm-Cm*Dm{3});
    >>> disp(err);
        4.656e-15
    >>> H = MinimalRepFromMRAP(Dm, 'obs');
    >>> disp(H{1});
          -4.4074       1.6931
          0.84259      -2.5926
    >>> disp(H{2});
          0.40741      0.13545
          0.55741     -0.20741
    >>> disp(H{3});
           1.6296       0.5418
           2.2296     -0.82963
    >>> Cm = SimilarityMatrix(H{1}, Dm{1});
    >>> err = norm(H{1}*Cm-Cm*Dm{1})+norm(H{2}*Cm-Cm*Dm{2})+norm(H{3}*Cm-Cm*Dm{3});
    >>> disp(err);
       4.8469e-15
    >>> H = MinimalRepFromMRAP(Dm, 'obscont');
    >>> disp(H{1});
          -4.4074       1.6931
          0.84259      -2.5926
    >>> disp(H{2});
          0.40741      0.13545
          0.55741     -0.20741
    >>> disp(H{3});
           1.6296       0.5418
           2.2296     -0.82963
    >>> Cm = SimilarityMatrix(H{1}, Dm{1});
    >>> err = norm(H{1}*Cm-Cm*Dm{1})+norm(H{2}*Cm-Cm*Dm{2})+norm(H{3}*Cm-Cm*Dm{3});
    >>> disp(err);
       4.8469e-15

    For Mathematica:

    
    For Python/Numpy:

    >>> D0 = ml.matrix([[-5., 1., 0],[3., -3., 0],[1., 1., -5.]])
    >>> D1 = ml.matrix([[0, 0, 0.8],[0, 0, 0],[0.2, 0.2, 0.2]])
    >>> D2 = ml.matrix([[0, 0, 3.2],[0, 0, 0],[0.8, 0.8, 0.8]])
    >>> Dm = [D0, D1, D2]
    >>> H = MinimalRepFromMRAP(Dm, "cont")
    >>> print(H[0])
    [[-5.  1.  0.]
     [ 3. -3.  0.]
     [ 1.  1. -5.]]
    >>> print(H[1])
    [[ 0.   0.   0.8]
     [ 0.   0.   0. ]
     [ 0.2  0.2  0.2]]
    >>> print(H[2])
    [[ 0.   0.   3.2]
     [ 0.   0.   0. ]
     [ 0.8  0.8  0.8]]
    >>> Cm = SimilarityMatrix(H[0], Dm[0])
    >>> err = la.norm(H[0]*Cm-Cm*Dm[0])+la.norm(H[1]*Cm-Cm*Dm[1])+la.norm(H[2]*Cm-Cm*Dm[2])
    >>> print(err)
    3.65128686018e-15
    >>> H = MinimalRepFromMRAP(Dm, "obs")
    >>> print(H[0])
    [[-4.40741  1.69312]
     [ 0.84259 -2.59259]]
    >>> print(H[1])
    [[ 0.40741  0.13545]
     [ 0.55741 -0.20741]]
    >>> print(H[2])
    [[ 1.62963  0.5418 ]
     [ 2.22963 -0.82963]]
    >>> Cm = SimilarityMatrix(H[0], Dm[0])
    >>> err = la.norm(H[0]*Cm-Cm*Dm[0])+la.norm(H[1]*Cm-Cm*Dm[1])+la.norm(H[2]*Cm-Cm*Dm[2])
    >>> print(err)
    4.72728076906e-15
    >>> H = MinimalRepFromMRAP(Dm, "obscont")
    >>> print(H[0])
    [[-4.40741  1.69312]
     [ 0.84259 -2.59259]]
    >>> print(H[1])
    [[ 0.40741  0.13545]
     [ 0.55741 -0.20741]]
    >>> print(H[2])
    [[ 1.62963  0.5418 ]
     [ 2.22963 -0.82963]]
    >>> Cm = SimilarityMatrix(H[0], Dm[0])
    >>> err = la.norm(H[0]*Cm-Cm*Dm[0])+la.norm(H[1]*Cm-Cm*Dm[1])+la.norm(H[2]*Cm-Cm*Dm[2])
    >>> print(err)
    4.72728076906e-15

