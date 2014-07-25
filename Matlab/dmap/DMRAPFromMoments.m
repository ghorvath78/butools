%  H = DMRAPFromMoments(moms, Nm)
%  
%  Creates a discrete marked rational arrival process that
%  has the same marginal and lag-1 joint moments as given 
%  (see [1]_).
%  
%  Parameters
%  ----------
%  moms : vector of doubles
%      The list of marginal moments. To obtain a discrete 
%      marked rational process of order M, 2*M-1 marginal 
%      moments are required.
%  Nm : list of matrices, shape (M,M)
%      The list of lag-1 joint moment matrices. The 
%      length of the list determines K, the number of arrival 
%      types of the discrete rational process.
%  
%  Returns
%  -------
%  H : list of matrices, shape (M,M)
%      The H0, H1, ..., HK matrices of the discrete marked
%      rational process
%  
%  References
%  ----------
%  .. [1] Andras Horvath, Gabor Horvath, Miklos Telek, "A 
%         traffic based decomposition of two-class queueing
%         networks with priority service," Computer Networks 
%         53:(8) pp. 1235-1248. (2009)

function H = DMRAPFromMoments(moms, Nm)

    [v, H0] = MGFromMoments (moms);
    N = size(H0,1);

    H0i = inv(eye(N)-H0);
    Ge = zeros(N);
    G1 = zeros(N);

    H0ip = eye(N);
    for i=1:N
        Ge(i,:) = v * H0ip;
        G1(:,i) = sum(H0ip, 2);
        H0ip = H0ip * i * H0i;
        if i>1
            H0ip = H0ip * H0;
        end
    end

    Gei = inv(Ge);
    G1i = inv(G1);   
    
    H = cell(1,length(Nm)+1);
    H{1} = H0;
    for i=2:length(Nm)+1
        Nmi = Nm{i-1};        
        row1 = FactorialMomsFromMoms(Nmi(1,2:end));
        col1 = FactorialMomsFromMoms(Nmi(2:end,1));
        mid = JFactorialMomsFromJMoms(Nmi(2:end,2:end));
        Nmi = [Nmi(1,1), row1; col1, mid];       
        H{i} = (eye(N)-H0)*Gei*Nmi*G1i;
    end
end
