%  Nm = LagkJointMomentsFromRAP(H0, H1, K, L, prec)
%  
%  Returns the lag-L joint moments of a rational arrival
%  process.
%  
%  Parameters
%  ----------
%  H0 : matrix, shape (M,M)
%      The H0 matrix of the rational arrival process
%  H1 : matrix, shape (M,M)
%      The H1 matrix of the rational arrival process
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
%  Nm : matrix, shape(K+1,K+1)
%      The matrix containing the lag-L joint moments, 
%      starting from moment 0.

function Nm = LagkJointMomentsFromRAP (H0, H1, K, L)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckRAPRepresentation(H0,H1)
        error('LagkJointMomentsFromRAP: Input isn''t a valid RAP representation!');
    end

    if ~exist('K','var')
        K=size(H0,1)-1;
    end

    if ~exist('L','var')
        L=1;
    end

    mom=LagkJointMomentsFromMRAP({H0, H1},K,L);
    Nm=mom{1};
end
