%  [H0, H1] = RAPFromMomentsAndCorrelations(moms, corr)
%  
%  Returns a rational arrival process that has the same moments
%  and lag autocorrelation coefficients as given.
%  
%  Parameters
%  ----------
%  moms : vector of doubles
%      The vector of marginal moments. To obtain a RAP of 
%      size M, 2*M-1 moments are required.
%  corr : vector of doubles
%      The vector of lag autocorrelation coefficients. To 
%      obtain a RAP of size M, 2*M-3 coefficients are needed.
%  
%  Returns
%  -------
%  H0 : matrix, shape (M,M)
%      The H0 matrix of the rational arrival process
%  H1 : matrix, shape (M,M)
%      The H1 matrix of the rational arrival process
%  
%  Notes
%  -----
%  There is no guarantee that the returned matrices define
%  a valid stochastic process. The joint densities may be
%  negative.
%  
%  References
%  ----------
%  .. [1] Mitchell, Kenneth, and Appie van de Liefvoort. 
%         "Approximation models of feed-forward G/G/1/N 
%         queueing networks with correlated arrivals." 
%         Performance Evaluation 51.2 (2003): 137-152.

function [D0,D1] = RAPFromMomentsAndCorrelations (moms, corr)

    [alpha, D0] = MEFromMoments (moms);
    M = length(alpha);
    
    if length(corr) < 2*M-3
        error('RAPFromMomentsAndCorrelations: The number of correlations given is less than required the 2n-3!');
    end
    
    rcorr=corr(1:2*M-3)/((moms(2)/2-moms(1)^2)/(moms(2)-moms(1)^2));
    rcorr = MomsFromReducedMoms (rcorr);
    [~,X] = MEFromMoments (reshape(rcorr,1,length(rcorr)));
    
    N = size(X,1);

    if N+1 ~= size(D0,1)
        error('RAPFromMomentsAndCorrelations: Correlation order is different from ME order');
    end

    T1 = zeros(N);
    for i=1:N
        for j=1:i
            T1(i,j) = 1;
        end
    end

    U1 = zeros(N);
    for i=1:N
        for j=i:N
            U1(i,j) = 1 / (N-i+1);
        end
    end

    T2 = zeros(M);
    for i=1:M
        for j=1:i
            T2(i,j) = 1;
        end
    end

    U2 = zeros(M);
    for i=1:M
        for j=i:M
            U2(i,j) = 1 / (M-i+1);
        end
    end

    Y = -inv(T1)*U1*inv(X)*inv(U1)*T1;
    II = eye(N+1);
    II(2:end,2:end)=Y;
    Y=II;

    D1=-D0*inv(U2)*T2*Y*inv(T2)*U2;
end
        





