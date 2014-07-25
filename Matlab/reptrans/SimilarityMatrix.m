%  B = SimilarityMatrix(A1, A2)
%  
%  Returns the matrix that transforms A1 to A2.
%  
%  Parameters
%  ----------
%  A1 : matrix, shape (N,N)
%      The smaller matrix
%  A2 : matrix, shape (M,M)
%      The larger matrix (M>=N)
%  
%  Returns
%  -------
%  B : matrix, shape (N,M)
%      The matrix satisfying `A_1\,B = B\,A_2`
%      
%  Notes
%  -----
%  For the existence of a (unique) solution the larger 
%  matrix has to inherit the eigenvalues of the smaller one.

function B = SimilarityMatrix (A1, A2)

    if size(A1,1)~=size(A1,2) || size(A2,1)~=size(A2,2)
        error('SimilarityMatrix: The input matrices must be square!');
    end

    if size(A1,1)>size(A2,1)
        error('SimilarityMatrix: The first input matrix must be smaller than the second one!');
    end

    N = size(A1,2);
    M = size(A2,1);
    X = kron(eye(M), A1) - kron(A2',eye(N));
    Y = kron(ones(1,M),eye(N));
    x = zeros(N*M,1);
    y = ones(N,1);
    Z = [X; Y];
    z = [x; y];
    B = reshape(linsolve(Z,z), N, M);
end
