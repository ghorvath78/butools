%  Nm = LagkJointMomentsFromMRAP(H, K, L, prec)
%  
%  Returns the lag-L joint moments of a marked rational 
%  arrival process.
%  
%  Parameters
%  ----------
%  H : list/cell of matrices of shape(M,M), length(N)
%      The H0...HN matrices of the MRAP to check
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
%  Nm : list/cell of matrices of shape(K+1,K+1), length(N)
%      The matrices containing the lag-L joint moments,
%      starting from moment 0.

function Nm = LagkJointMomentsFromMRAP (H, K, L)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMRAPRepresentation(H)
        error('LagkJointMomentsFromMRAP: Input isn''t a valid MRAP representation!');
    end

    M=size(H,2)-1;

    if ~exist('K','var')
        K=size(H{1},1)-1;
    end

    if ~exist('L','var')
        L=1;
    end

    sumH = SumMatrixList(H(2:end));
    iH0= inv(-H{1});

    pi = DRPSolve(iH0*sumH);
    Pw = eye(size(H{1},1));
    H0p=cell(K+1);
    H0p{1} = Pw;
    for i=1:K
        Pw = Pw*i*iH0;
        H0p{i+1} = Pw;
    end

    Pl = (iH0*sumH)^(L-1);
    Nm=cell(1,M);

    for m=1:M
        Nmm = zeros(K+1,K+1);
        for i=1:K+1
            for j=1:K+1
                Nmm(i,j)=sum(pi*H0p{i}*iH0*H{m+1}*Pl*H0p{j});
            end
        end
        Nm{m}=Nmm;
    end
end
