%  Nm = LagkJointMomentsFromDMRAP(H, K, L, prec)
%  
%  Returns the lag-L joint moments of a discrete marked 
%  rational arrival process.
%  
%  Parameters
%  ----------
%  H : list/cell of matrices of shape(M,M), length(N)
%      The H0...HN matrices of the DMRAP to check
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

function [Nm] = LagkJointMomentsFromDMRAP (h, K, L)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDMRAPRepresentation(h)
        error('LagkJointMomentsFromDMRAP: Input isn''t a valid DMRAP representation!');
    end

    M=size(h,2)-1;

    if ~exist('K','var')
        K=size(h{1},1)-1;
    end

    if ~exist('L','var')
        L=1;
    end

    sumH = SumMatrixList(h(2:end));
    iH0= inv(eye(size(h{1},1))-h{1});

    pi = DRPSolve(iH0*sumH);
    Pw = eye(size(h{1},1));
    H0p = cell(K+1);
    H0p{1} = Pw;
    Pw = Pw * iH0;
    H0p{2} = Pw;
    for i=2:K
        Pw = Pw*i*iH0*h{1};
        H0p{i+1} = Pw;
    end

    Pl = (iH0*sumH)^(L-1);
    Nm=cell(1,M);

    for m=1:M
        Nmm = zeros(K+1,K+1);
        for i=1:K+1
            for j=1:K+1
                Nmm(i,j)=sum(pi*H0p{i}*iH0*h{m+1}*Pl*H0p{j});
            end
        end
        row1 = MomsFromFactorialMoms(Nmm(1,2:end));
        col1 = MomsFromFactorialMoms(Nmm(2:end,1));
        mid = JMomsFromJFactorialMoms(Nmm(2:end,2:end));
        Nm{m}= [Nmm(1,1), row1; col1, mid];
    end

end