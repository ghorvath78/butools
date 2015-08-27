%  [mass0, ini, K, clo] = FluidSolve (Fpp, Fpm, Fmp, Fmm, prec)
%  
%  Returns the parameters of the matrix-exponentially 
%  distributed stationary distribution of a canonical 
%  Markovian fluid model.
%  
%  The canonical Markov fluid model is defined by the 
%  matrix blocks of the generator of the background Markov
%  chain partitioned according to the sign of the 
%  associated fluid rates (i.e., there are "+" and "-" states).   
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

function [mass0, ini, K, clo] = FluidSolve (Fpp, Fpm, Fmp, Fmm, prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end

    [Psi, K, U] = FluidFundamentalMatrices (Fpp, Fpm, Fmp, Fmm, 'PKU', prec);
    mass0 = CTMCSolve(U);
    nr = sum(mass0) + 2*sum(mass0*Fmp*inv(-K));
    mass0 = mass0 / nr;
    ini = mass0 * Fmp;
    clo = [eye(size(Fpp,1)), Psi];
    mass0 = [zeros(1,size(Fpp,1)), mass0];       
end

