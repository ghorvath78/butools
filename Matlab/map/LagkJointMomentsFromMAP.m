%  Nm = LagkJointMomentsFromMAP(D0, D1, K, L, prec)
%  
%  Returns the lag-L joint moments of a Markovian arrival
%  process.
%  
%  Parameters
%  ----------
%  D0 : matrix, shape (M,M)
%      The D0 matrix of the Markovian arrival process
%  D1 : matrix, shape (M,M)
%      The D1 matrix of the Markovian arrival process
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

function Nm = LagkJointMomentsFromMAP (D0, D1, K, L)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMAPRepresentation(D0,D1)
        error('LagkJointMomentsFromMAP: Input isn''t a valid MAP representation!');
    end

    if ~exist('K','var')
        K=size(D0,1)-1;
    end

    if ~exist('L','var')
        L=1;
    end

    moms=LagkJointMomentsFromMRAP({D0,D1},K,L);
    Nm=moms{1};
end
