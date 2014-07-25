butools.queues.FluidQueueSTD
============================

.. currentmodule:: butools.queues

.. np:function:: FluidQueueSTD

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = FluidQueueSTD (Q, Rin, Rout, Q0, transToPH, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = FluidQueueSTD [Q, Rin, Rout, Q0, transToPH, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = FluidQueueSTD (Q, Rin, Rout, Q0, transToPH, prec)`

    Returns the parameters of the matrix-exponentially 
    distributed stationary sojourn time distribution of
    fluid drops in a fluid queue.
    
    In a fluid queue there is a background continuous time
    Markov chain (given by generator Q), and diagonal
    matrix Rin (Rout) whose ith entry provides the 
    fluid rate at which fluid enters the queue (can be 
    served) while the background process is in state i.
    
    
    Parameters
    ----------
    Q : matrix, shape (N,N)
        The generator of the background Markov chain
    Rin : matrix, shape (N,N)
        Diagonal matrix containing the fluid input rates
        associated to the states of the background process
    Rout : matrix, shape (N,N)
        Diagonal matrix containing the fluid output rates
        associated to the states of the background process
    Q0 : matrix, shape (N,N), optional
        The generator of the background Markov chain when
        the queue level is zero. If empty or missing, Q0=Q
        is assumed
    transToPH : bool, optional
        If true, the result is transformed to a phase-type
        representation, otherwise it is returned as a matrix
        exponential representation. The default value is 
        true
    prec : double, optional
        Numerical precision used to check wether the input 
        is valid and it also serves as a stopping condition
        when the algebraic Riccati equation is solved. The
        delault value is 1e-14
    
    Returns
    -------
    alpha : matrix, shape (1,M)
        The initial vector of the matrix-exponential 
        sojourn time distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential 
        sojourn time distribution.
    
    Notes
    -----
    While it is always possible to transform the result to 
    a phase-type representation theoretically, this 
    transformation step can be sensitive numerically.
    
    Notes
    -----
    The result can be very large! If the input and output
    processes are independent, use the 'FluFluSTD' function
    which much faster and returns a small representation.

    Examples
    ========
    For MATLAB:

    >>> Q = [-9 2 4 0 1 2; 6 -25 5 3 7 4; 1 3 -4 0 0 0; 0 0 0 -8 3 5; 7 3 0 2 -13 1; 7 8 0 3 8 -26];
    >>> Rin = diag([4 2 1 0 0 3]);
    >>> Rout = diag([6 2 0 0 3 2]);
    >>> [alpha, A] = FluidQueueQLD(Q,Rin,Rout);
    >>> MomentsFromME(alpha,A,5)
          0.23252      0.20069      0.26684      0.47402       1.0523
    >>> x=(0:0.05:1)';
    >>> y=PdfFromME(alpha,A,x);
    >>> plot(x,y);
    >>> [alpha, A] = FluidQueueSTD(Q,Rin,Rout,[],true);
    >>> CheckPHRepresentation(alpha,A)
         1    
            
