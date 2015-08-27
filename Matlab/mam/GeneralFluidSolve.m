%  [mass0, ini, K, clo] = GeneralFluidSolve (Q, R, Q0, prec)
%  
%  Returns the parameters of the matrix-exponentially 
%  distributed stationary distribution of a general 
%  Markovian fluid model, where the fluid rates associated
%  with the states of the background process can be
%  arbitrary (zero is allowed as well).
%  
%  Using the returned 4 parameters the stationary
%  solution can be obtained as follows.
%  
%  The probability that the fluid level is zero while 
%  being in different states of the background process
%  is given by vector mass0.
%  
%  The density that the fluid level is x while being in
%  different states of the background process is
%  
%  .. math::
%      \pi(x)=ini\cdot e^{K x}\cdot clo.    
%  
%  Parameters
%  ----------
%  Q : matrix, shape (N,N)
%      The generator of the background Markov chain
%  R : diagonal matrix, shape (N,N)
%      The diagonal matrix of the fluid rates associated
%      with the different states of the background process
%  Q0 : matrix, shape (N,N), optional
%      The generator of the background Markov chain at 
%      level 0. If not provided, or empty, then Q0=Q is 
%      assumed. The default value is empty.
%  precision : double, optional
%      Numerical precision for computing the fundamental
%      matrix. The default value is 1e-14
%  
%  Returns
%  -------
%  mass0 : matrix, shape (1,Np+Nm)
%      The stationary probability vector of zero level
%  ini : matrix, shape (1,Np)
%      The initial vector of the stationary density
%  K : matrix, shape (Np,Np)
%      The matrix parameter of the stationary density
%  clo : matrix, shape (Np,Np+Nm)
%      The closing matrix of the stationary density

function [mass0, ini, K, clo] = GeneralFluidSolve (Q, R, Q0, prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end

    if ~exist('Q0','var')
        Q0 = [];
    end
    
    N = size(Q,1);

    % partition the state space according to zero, positive and negative fluid rates
    ix = (1:N);
    ixz = ix(abs(diag(R))<=prec);
    ixp = ix(diag(R)>prec);
    ixn = ix(diag(R)<-prec);
    Nz = length(ixz);
    Np = length(ixp);
    Nn = length(ixn);
    
    % permutation matrix that converts between the original and the partitioned state ordering
    P = zeros(N);
    for i=1:Nz
        P(i,ixz(i))=1;
    end
    for i=1:Np
        P(Nz+i,ixp(i))=1;
    end
    for i=1:Nn
        P(Nz+Np+i,ixn(i))=1;
    end
    iP = inv(P);
    
    % reorder states with permutation matrix P: 0, +, -
    Qv = P*Q*iP;
    Rv = P*R*iP;
    
    % new fluid process censored to states + and -
    iQv00 = pinv(-Qv(1:Nz,1:Nz));
    Qbar = Qv(Nz+1:end, Nz+1:end) + Qv(Nz+1:end,1:Nz)*iQv00*Qv(1:Nz,Nz+1:end);   
    absRi = diag(abs(1./diag(Rv(Nz+1:end,Nz+1:end))));
    Qz = absRi * Qbar;
    
    % calculate fundamental matrices
    [Psi, K, U] = FluidFundamentalMatrices (Qz(1:Np,1:Np), Qz(1:Np,Np+1:end), Qz(Np+1:end,1:Np), Qz(Np+1:end,Np+1:end), 'PKU', prec);
    
    % closing matrix
    Pm = [eye(Np), Psi];
    iCn = absRi(Np+1:end,Np+1:end);
    iCp = absRi(1:Np,1:Np);
    clo = [(iCp*Qv(Nz+1:Nz+Np,1:Nz)+Psi*iCn*Qv(Nz+Np+1:end,1:Nz))*iQv00, Pm*absRi];

    if isempty(Q0) % regular boundary behavior
    
        clo = clo * P; % go back the the original state ordering
        
        % calculate boundary vector   
        Ua = iCn*Qv(Nz+Np+1:end,1:Nz)*iQv00*ones(Nz,1) + iCn*ones(Nn,1) + Qz(Np+1:end,1:Np)*inv(-K)*clo*ones(Nz+Np+Nn,1);
        pm = linsolve ([U,Ua]', [zeros(1,Nn),1]')';

        % create the result
        mass0 = [pm*iCn*Qv(Nz+Np+1:end,1:Nz)*iQv00, zeros(1,Np), pm*iCn]*P;
        ini = pm*Qz(Np+1:end,1:Np);
    else
        
        % solve a linear system for ini(+), pm(-) and pm(0)        
        Q0v = P*Q0*iP;
        M = [-clo*Rv; Q0v(Nz+Np+1:end,:); Q0v(1:Nz,:)];
        Ma = [sum(inv(-K)*clo,2); ones(Nz+Nn,1)];
        sol = linsolve ([M,Ma]', [zeros(1,N),1]')';
        ini = sol(1:Np);
        clo = clo * P;
        mass0 = [sol(Np+Nn+1:end), zeros(1,Np), sol(Np+1:Np+Nn)]*P;
    end
end