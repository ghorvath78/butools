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

    N1 = size(A1,1);
    N2 = size(A2,1);

    if N1>N2
        error('SimilarityMatrix: The first input matrix must be smaller than the second one!');
    end
       
    Ax = A2;
    Ax(:,1)=Ax(:,1)+ones(N2,1);
    C = zeros(N1,N2);
    C(:,1) = C(:,1)+ones(N1,1);
    B = lyap(A1,-Ax,C);    
end
