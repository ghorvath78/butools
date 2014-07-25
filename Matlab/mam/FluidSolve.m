%  [mass0, ini, K, clo] = FluidSolve (Fpp, Fpm, Fmp, Fmm)
%  
%  Returns the parameters of the matrix-exponentially 
%  distributed stationary distribution of a canonical 
%  Markovian fluid model
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
%  Fpp : matrix, shape (Np,Np)
%      The matrix of transition rates between states 
%      having positive fluid rates
%  Fpm : matrix, shape (Np,Nm)
%      The matrix of transition rates where the source
%      state has a positive, the destination has a 
%      negative fluid rate associated.
%  Fpm : matrix, shape (Nm,Np)
%      The matrix of transition rates where the source
%      state has a negative, the destination has a 
%      positive fluid rate associated.
%  Fpp : matrix, shape (Nm,Nm)
%      The matrix of transition rates between states 
%      having negative fluid rates
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

function [mass0, ini, K, clo] = FluidSolve (Fpp, Fpm, Fmp, Fmm, prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end

    M = FluidFundamentalMatrices (Fpp, Fpm, Fmp, Fmm, 'PKU', prec);
    Psi = M{1};
    K = M{2};
    U = M{3};    
    mass0 = CTMCSolve(U, prec);
    nr = sum(mass0) + 2*sum(mass0*Fmp*inv(-K));
    mass0 = mass0 / nr;       
    ini = mass0 * Fmp;
    clo = [eye(size(Fpp,1)), Psi];
end

