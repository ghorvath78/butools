%  [alpha, A] = FluFluSTD (Qin, Rin, Qout, Rout, srv0stop, transToPH, prec)
%  
%  Returns the parameters of the matrix-exponentially 
%  distributed stationary sojourn time distribution of
%  fluid drops in a fluid queue where the fluid input
%  and output processes are independent.
%  
%  Two types of boundary behavior is available. If 
%  srv0stop=false, the output process evolves continuously
%  even if the queue is empty. If srv0stop=true, the 
%  output process slows down if there is fewer fluid in
%  the queue than it can serve. If the queue is empty
%  and the fluid input rate is zero, the output process
%  freezes till fluid arrives.
%  
%  Parameters
%  ----------
%  Qin : matrix, shape (N,N)
%      The generator of the background Markov chain 
%      corresponding to the input process
%  Rin : matrix, shape (N,N)
%      Diagonal matrix containing the fluid input rates
%      associated to the states of the input background 
%      process
%  Qout : matrix, shape (N,N)
%      The generator of the background Markov chain 
%      corresponding to the output process
%  Rout : matrix, shape (N,N)
%      Diagonal matrix containing the fluid output rates
%      associated to the states of the input background 
%      process
%  transToPH : bool, optional
%      If true, the result is transformed to a phase-type
%      representation, otherwise it is returned as a matrix
%      exponential representation. The default value is 
%      true
%  prec : double, optional
%      Numerical precision used to check wether the input 
%      is valid and it also serves as a stopping condition
%      when the algebraic Riccati equation is solved. The
%      delault value is 1e-14
%  
%  Returns
%  -------
%  alpha : matrix, shape (1,M)
%      The initial vector of the matrix-exponential 
%      sojourn time distribution.
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix-exponential 
%      sojourn time distribution.
%  
%  Notes
%  -----
%  While it is always possible to transform the result to 
%  a phase-type representation theoretically, this 
%  transformation step can be sensitive numerically.

function [alpha, A] = FluFluSTD (Qin, Rin, Qout, Rout, srv0stop, transToPH, prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end

    if ~exist('srv0stop','var')
        srv0stop = false;
    end

    if ~exist('transToPH','var')
        transToPH = false;
    end

    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckGenerator(Qin,false,prec)
        error('FluFluQLD: Generator matrix Qin is not Markovian!');
    end
    if BuToolsCheckInput && ~CheckGenerator(Qout,false,prec)
        error('FluFluQLD: Generator matrix Qout is not Markovian!');
    end
    if BuToolsCheckInput && (any(diag(Rin)<-prec) || any(diag(Rout)<-prec))
        error('FluFluQLD: Fluid rates Rin and Rout must be non-negative !');
    end
    
    % solve special fluid queue
    Iin = eye(size(Qin));
    Iout = eye(size(Qout));

    Rh = kron(Rin,Iout) - kron(Iin,Rout);
    Qh = kron(Qin, Rout) + kron(Rin, Qout);       
    [massh, inih, Kh, cloh] = GeneralFluidSolve (Qh, Rh);

    % sojourn time density in case of 
    % srv0stop = false: inih*expm(Kh*x)*cloh*kron(Rin,Iout)/lambda
    % srv0stop = true: inih*expm(Kh*x)*cloh*kron(Rin,Rout)/lambda/mu    
    
    lambda = sum(CTMCSolve(Qin)*Rin);
    mu = sum(CTMCSolve(Qout)*Rout);
    
    if transToPH   
        % convert result to PH representation
        Delta = diag(linsolve(Kh',-inih')); % Delta = diag (inih*inv(-Kh));
        A = inv(Delta)*Kh'*Delta;       
        if ~srv0stop        
            alpha = sum(Delta*cloh*kron(Rin,Iout)/lambda,2)';
        else
            alpha = sum(Delta*cloh*kron(Rin,Rout)/lambda/mu,2)';
        end        
    else
        % convert result to ME representation
        if ~srv0stop
            B = TransformToOnes(sum(cloh*kron(Rin,Iout)/lambda,2));
        else
            B = TransformToOnes(sum(cloh*kron(Rin,Rout)/lambda/mu,2));
        end
        iB = inv(B);
        A = B*Kh*iB;
        alpha = inih*inv(-Kh)*iB;
    end
end

