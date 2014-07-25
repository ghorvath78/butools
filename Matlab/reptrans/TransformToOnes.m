%  B = TransformToOnes(clovec)
%  
%  Returns the similarity transformation matrix that converts 
%  the given column vector to a vector of ones. It works even
%  if it has zero entries.
%  
%  Parameters
%  ----------
%  clovec : column vector, shape(M,1)
%      The original closing vector
%      
%  Returns
%  -------
%  B : matrix, shape(M,M)
%      The matrix by which `B\cdot clovec = \mathbf{1}` holds

function B = TransformToOnes (clovec)

    % construct permutation matrix to move at least one non-zero element to the first position
    % to acchieve it, the code below sorts it in a reverse order
    m = length(clovec);
    
    [~,ix] = sort(-clovec);
    P = zeros(m);
    for i=1:m
        P(i,ix(i)) = 1.0;
    end
    cp = P*clovec;

    % construct transformation matrix B for which B*rp=1 holds
    B = zeros(m);
    for i=1:m
        B(i,1:i) = 1.0 / sum(cp(1:i,1));
    end
    % construct matrix B for which B*r=1 holds
    B = B*P;
end
