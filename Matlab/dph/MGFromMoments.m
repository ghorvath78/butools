%  [alpha, A] = MGFromMoments(moms)
%  
%  Creates a matrix-geometric distribution that has the
%  same moments as given.
%  
%  Parameters
%  ----------
%  moms : vector of doubles
%      The list of moments. The order of the resulting 
%      matrix-geometric distribution is 
%      determined based on the number of moments given. To 
%      obtain a matrix-geometric distribution of order M,
%      2*M-1 moments are required.
%  
%  Returns
%  -------
%  alpha : vector, shape (1,M)
%      The initial vector of the matrix-geometric 
%      distribution.
%  A : matrix, shape (M,M)
%      The matrix parameter of the matrix-geometric 
%      distribution.
%  
%  References
%  ----------
%  .. [1] A. van de Liefvoort. The moment problem for 
%         continuous distributions. Technical report, 
%         University of Missouri, WP-CM-1990-02, Kansas City,
%         1990.

function [alpha,A] = MGFromMoments(moms)

    rfmoms = ReducedMomsFromMoms (FactorialMomsFromMoms(moms));
    rfmoms = [1,rfmoms];

    vlist = zeros(1,length(moms));
    tmpVec=zeros(1,length(moms)+1);
    k=1;
    tmpVec(1)=rfmoms(1);

    for i=1:length(moms)
        tmpVec(i+1)=(-1)^i*rfmoms(i+1);
        k=k*i;
        vlist(i)=k*sum(tmpVec);
    end

    [alpha,C]=MEFromMoments(vlist);
    iC=inv(C);
    A=iC*inv(iC+eye(length(iC)));
end
