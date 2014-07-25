%  H = MRAPFromMoments(moms, Nm)
%  
%  Creates a marked rational arrival process that has the same 
%  marginal and lag-1 joint moments as given (see [1]_).
%  
%  Parameters
%  ----------
%  moms : vector of doubles
%      The list of marginal moments. To obtain a marked 
%      rational process of order M, 2*M-1 marginal moments
%      are required.
%  Nm : list of matrices, shape (M,M)
%      The list of lag-1 joint moment matrices. The 
%      length of the list determines K, the number of arrival 
%      types of the rational process.
%  
%  Returns
%  -------
%  H : list of matrices, shape (M,M)
%      The H0, H1, ..., HK matrices of the marked rational 
%      process
%  
%  Notes
%  -----
%  There is no guarantee that the returned matrices define
%  a valid stochastic process. The joint densities may be
%  negative.
%  
%  References
%  ----------
%  .. [1] Andras Horvath, Gabor Horvath, Miklos Telek, "A 
%         traffic based decomposition of two-class queueing
%         networks with priority service," Computer Networks 
%         53:(8) pp. 1235-1248. (2009)

function H = MRAPFromMoments (moms, Nm)

    [v, H0] = MEFromMoments (moms);
    H0i = inv(-H0);

    N = size(H0,1);
    Ge = zeros(N);
    G1 = zeros(N);

    H0ip = eye(N);
    for i=1:N
        Ge(i,:) = v * H0ip;
        G1(:,i) = sum(H0ip, 2);
        H0ip = H0ip * i * H0i;
    end

    Gei = inv(Ge);
    G1i = inv(G1);
    
    H = cell(1,length(Nm)+1);
    H{1} = H0;
    for i=2:length(Nm)+1
        H{i} = -H0*Gei*Nm{i-1}*G1i;
    end
end

