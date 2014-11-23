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
    
    [Q1,R1]=schur(A1,'complex');
    [Q2,R2]=schur(A2,'complex');
    
    c1 = sum(Q2',2);
    c2 = sum(Q1',2);
    I = eye(N2);
    X = zeros(N1,N2);
    for k=N1:-1:1
        M = R1(k,k)*I-R2;
        if k==N1
            m = zeros(1,N2);
        else
            m = -R1(k,k+1:end)*X(k+1:end,:);
        end
        X(k,:) = linsolve([M,c1]',[m,c2(k)]')';
    end    
    B = real(Q1*X*Q2');
end
