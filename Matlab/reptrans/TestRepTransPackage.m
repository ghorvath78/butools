format short g;

disp('---BuTools: RepTrans package test file---');

disp('Enable the verbose messages with the BuToolsVerbose flag')
global BuToolsVerbose;
BuToolsVerbose = true

disp('Enable input parameter checking with the BuToolsCheckInput flag')
global BuToolsCheckInput;
BuToolsCheckInput = true

disp('----------------------------------------------------------------------------');
help SimilarityMatrix

disp('Input:');
disp('------');

A1 = [0.2 0.8 0; 1.2 -0.4 0.1; -0.2 0.7 0.5]
T = [1 2 -4 6; 0 8 -9 7; -3 7 8 -2];
A2 = pinv(T)*A1*T
%A2 = [0.0088058 0.29904 -0.4672 0.24887; -0.034691 -0.0035465 0.45729 0.02035; -0.20158 1.7687 -1.3378 0.82115; -0.090959 2.2638 -2.2997 1.6325]

disp('Test:');
disp('-----');

disp('B=SimilarityMatrix(A1,A2):');
B=SimilarityMatrix(A1,A2);
disp(B);
err = norm(A1*B-B*A2)
assert(err<sqrt(eps), 'The resulting matrix T does not satisfy A1*T = T*A2!');

disp('----------------------------------------------------------------------------');
help TransformToAcyclic

disp('Input:');
disp('------');

A = [-0.8 0.8 0; 0.1 -0.3 0.1; 0.2 0 -0.5]

disp('Test:');
disp('-----');

disp('B=TransformToAcyclic(A):');
B=TransformToAcyclic(A);
disp(B);
C=SimilarityMatrix(A,B);
err = norm(A*C-C*B)
assert(err<sqrt(eps), 'The original and the transformed matrix are not similar!');

disp('----------------------------------------------------------------------------');
help TransformToMonocyclic

disp('Input:');
disp('------');

% A = [-0.7 0.7 0.0; 0.0 -0.8 0.8; 0.7 0.0 -0.8]
A = [-1,0,0;0,-3,2;0,-2,-3]

disp('Test:');
disp('-----');

disp('B=TransformToMonocyclic(A):');
B=TransformToMonocyclic(A);
disp(B);
C=SimilarityMatrix(A,B);
err = norm(A*C-C*B)
assert(err<sqrt(eps), 'The original and the transformed matrix are not similar!');

disp('----------------------------------------------------------------------------');
help ExtendToMarkovian

disp('Input:');
disp('------');

disp('Original PH:');
%alpha = [0.2 0.0 0.8]
%A = [-0.7 0.7 0.0; 0.0 -0.8 0.8; 0.7 0.0 -0.8]
alpha = [0.2, 0.3, 0.5]
A = [-1,0,0;0,-3,0.6;0,-0.6,-3] 

disp('Transformed to Monocyclic:');
B=TransformToMonocyclic(A);
disp(B);
C=SimilarityMatrix(A,B);
beta = alpha*C;
disp(beta)

disp('Test:');
disp('-----');

disp('[m,M]=ExtendToMarkovian(beta,B):');
[m,M]=ExtendToMarkovian(beta,B);
disp(m);
disp(M);
C=SimilarityMatrix(B,M);
err = norm(B*C-C*M)
assert(err<sqrt(eps), 'The original and the transformed matrix are not similar!');
assert(min(m)>-1e-14, 'The initial vector is still not Markovian!');

disp('----------------------------------------------------------------------------');
help SimilarityMatrixForVectors

disp('Input:');
disp('------');

disp('Original column vectors:');
vecA = [0.0, 0.3, -1.5, 0.0]'
vecB = [1.0, 0.2, 0.0, 1.0]'

disp('Test:');
disp('-----');

disp('B=SimilarityMatrixForVectors (vecA, vecB):');
B = SimilarityMatrixForVectors (vecA, vecB);
disp(B);

assert(max(abs(B*vecA-vecB))<1e-14, 'The resulting matrix T does not satisfy T*vecA = vecB!');

disp('----------------------------------------------------------------------------');
help FindMarkovianRepresentation

disp('Tested in the PH and in the MAP package');

disp('----------------------------------------------------------------------------');
help MStaircase

disp('Tested in the PH and in the MAP package');

