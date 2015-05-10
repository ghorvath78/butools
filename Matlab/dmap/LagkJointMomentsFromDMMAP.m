%  Nm = LagkJointMomentsFromDMMAP(D, K, L, prec)
%  
%  Returns the lag-L joint moments of a discrete marked 
%  Markovian arrival process.
%  
%  Parameters
%  ----------
%  D : list/cell of matrices of shape(M,M), length(N)
%      The D0...DN matrices of the DMMAP to check
%  K : int, optional
%      The dimension of the matrix of joint moments to 
%      compute. If K=0, the MxM joint moments will be 
%      computed. The default value is 0
%  L : int, optional
%      The lag at which the joint moments are computed.
%      The default value is 1
%  prec : double, optional
%      Numerical precision to check if the input is valid. 
%      The default value is 1e-14
%  
%  Returns
%  -------
%  Nm : list/cell of matrices of shape(K+1,K+1), length(L)
%      The matrices containing the lag-L joint moments,
%      starting from moment 0.

function [Nm] = LagkJointMomentsFromDMMAP (h, K, L)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDMMAPRepresentation(h)
        error('LagkJointMomentsFromDMMAP: Input isn''t a valid DMMAP representation!');
    end

    if ~exist('K','var')
        K=size(h{1},1)-1;
    end

    if ~exist('L','var')
        L=1;
    end

    Nm=LagkJointMomentsFromDMRAP(h,K,L);
end