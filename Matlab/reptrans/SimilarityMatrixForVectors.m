%  B = SimilarityMatrixForVectors(vecA, vecB)
%  
%  Returns the similarity transformation matrix that converts 
%  a given column vector to an other column vector. It works 
%  even with zero entries.
%  
%  Parameters
%  ----------
%  vecA : column vector, shape(M,1)
%      The original column vector
%  vecB : column vector, shape(M,1)
%      The target column vector
%      
%  Returns
%  -------
%  B : matrix, shape(M,M)
%      The matrix by which `B\cdot vecA = vecB` holds

function B = SimilarityMatrixForVectors (vecA, vecB)

    % construct permutation matrix to move at least one non-zero element to the first position
    % to acchieve it, the code below sorts it in a reverse order
    m = length(vecA);
    
    [~,ix] = sort(-vecA);
    P = zeros(m);
    for i=1:m
        P(i,ix(i)) = 1.0;
    end
    cp = P*vecA;

    % construct transformation matrix B for which B*rp=1 holds
    B = zeros(m);
    for i=1:m
        B(i,1:i) = vecB(i) / sum(cp(1:i,1));
    end
    % construct matrix B for which B*r=1 holds
    B = B*P;
end
