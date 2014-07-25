butools.queues.FluFluQLD
========================

.. currentmodule:: butools.queues

.. np:function:: FluFluQLD

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`[alpha, A] = FluFluQLD (Qin, Rin, Qout, Rout, srv0stop, transToPH, prec)`
        * - Mathematica:
          - :code:`{alpha, A} = FluFluQLD [Qin, Rin, Qout, Rout, srv0stop, transToPH, prec]`
        * - Python/Numpy:
          - :code:`alpha, A = FluFluQLD (Qin, Rin, Qout, Rout, srv0stop, transToPH, prec)`

    Returns the parameters of the matrix-exponentially 
    distributed stationary fluid level distribution of
    fluid drops in a fluid queue where the fluid input
    and output processes are independent.
    
    Two types of boundary behavior is available. If 
    srv0stop=false, the output process evolves continuously
    even if the queue is empty. If srv0stop=true, the 
    output process slows down if there is fewer fluid in
    the queue than it can serve. If the queue is empty
    and the fluid input rate is zero, the output process
    freezes till fluid arrives.
    
    Parameters
    ----------
    Qin : matrix, shape (N,N)
        The generator of the background Markov chain 
        corresponding to the input process
    Rin : matrix, shape (N,N)
        Diagonal matrix containing the fluid input rates
        associated to the states of the input background 
        process
    Qout : matrix, shape (N,N)
        The generator of the background Markov chain 
        corresponding to the output process
    Rout : matrix, shape (N,N)
        Diagonal matrix containing the fluid output rates
        associated to the states of the input background 
        process
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
        fluid level distribution.
    A : matrix, shape (M,M)
        The matrix parameter of the matrix-exponential 
        fluid level distribution.
    
    Notes
    -----
    While it is always possible to transform the result to 
    a phase-type representation theoretically, this 
    transformation step can be sensitive numerically.
 
    Examples
    ========
    For MATLAB:

    >>> Qin = [-2 1 1; 2 -5 3; 4 0 -4];
    >>> Rin = diag([3 7 0]);
    >>> Qout = [-4 1 3; 6 -8 2; 3 7 -10];
    >>> Rout = diag([1 7 15]);
    >>> [alpha, A] = FluFluQLD(Qin,Rin,Qout,Rout,false);
    >>> MomentsFromME(alpha,A,5)
           0.5357       1.0765       3.4298        14.87       81.162
    >>> x=(0:0.1:2)';
    >>> y=PdfFromME(alpha,A,x);
    >>> plot(x,y);
    >>> [alpha, A] = FluFluQLD(Qin,Rin,Qout,Rout,false,true);
    >>> CheckPHRepresentation(alpha,A)
         1    
        
