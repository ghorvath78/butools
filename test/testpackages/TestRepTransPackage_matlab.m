function TestRepTransPackage ()
	disp('---BuTools: RepTrans package test file---');
	disp('Enable the verbose messages with the BuToolsVerbose flag');
	global BuToolsVerbose;
	BuToolsVerbose = true;
	disp('Enable input parameter checking with the BuToolsCheckInput flag');
	global BuToolsCheckInput;
	BuToolsCheckInput = true;
	global BuToolsCheckPrecision;
	disp('========================================')
	help SimilarityMatrix
	A1m = [0.2, 0.8, 0; 1.2, -0.4, 0.1; -0.2, 0.7, 0.5];
	disp('A1m = ');
	disp(A1m);
	T = [1, 2, -4, 6; 0, 8, -9, 7; -3, 7, 8, -2];
	disp('T = ');
	disp(T);
	disp('A2m = pinv(T)*A1m*T:');
	A2m = pinv(T)*A1m*T;
	disp('Test:');
	disp('-----');
	disp('B = SimilarityMatrix(A1m,A2m):');
	B = SimilarityMatrix(A1m,A2m);
	disp('err = norm(A1m*B-B*A2m):');
	err = norm(A1m*B-B*A2m);
	disp('err = ');
	disp(err);
	assert(err<10^-7 , 'The resulting matrix T does not satisfy A1m*T = T*A2m!');
	disp('========================================')
	help TransformToAcyclic
	disp('Input:');
	disp('------');
	A = [-0.8, 0.8, 0; 0.1, -0.3, 0.1; 0.2, 0, -0.5];
	disp('A = ');
	disp(A);
	disp('Test:');
	disp('-----');
	disp('B = TransformToAcyclic(A):');
	B = TransformToAcyclic(A);
	disp('B = ');
	disp(B);
	disp('Cm = SimilarityMatrix(A,B):');
	Cm = SimilarityMatrix(A,B);
	disp('err = norm(A*Cm-Cm*B):');
	err = norm(A*Cm-Cm*B);
	disp('err = ');
	disp(err);
	assert(err<10^-7 , 'The original and the transformed matrix are not similar!');
	disp('========================================')
	help TransformToMonocyclic
	A = [-1, 0, 0; 0, -3, 2; 0, -2, -3];
	disp('A = ');
	disp(A);
	disp('Test:');
	disp('-----');
	disp('B = TransformToMonocyclic(A):');
	B = TransformToMonocyclic(A);
	disp('B = ');
	disp(B);
	disp('Cm = SimilarityMatrix(A,B):');
	Cm = SimilarityMatrix(A,B);
	disp('err = norm(A*Cm-Cm*B):');
	err = norm(A*Cm-Cm*B);
	disp('err = ');
	disp(err);
	assert(err<10^-7 , 'The original and the transformed matrix are not similar!');
	disp('========================================')
	help ExtendToMarkovian
	disp('Input:');
	disp('------');
	alpha = [0.2, 0.3, 0.5];
	A = [-1, 0, 0; 0, -3, 0.6; 0, -0.6, -3];
	B = TransformToMonocyclic(A);
	disp('B = ');
	disp(B);
	Cm = SimilarityMatrix(A,B);
	beta = alpha*Cm;
	disp('beta = ');
	disp(beta);
	disp('Test:');
	disp('-----');
	disp('[m,M] = ExtendToMarkovian(beta,B):');
	[m,M] = ExtendToMarkovian(beta,B);
	disp('m = ');
	disp(m);
	disp('M = ');
	disp(M);
	disp('Cm = SimilarityMatrix(B,M):');
	Cm = SimilarityMatrix(B,M);
	disp('err = norm(B*Cm-Cm*M):');
	err = norm(B*Cm-Cm*M);
	disp('err = ');
	disp(err);
	assert(err<10^-7 , 'The original and the transformed matrix are not similar!');
	assert(min(m)>-10^-14 , 'The initial vector is still not Markovian!');
	disp('========================================')
	help SimilarityMatrixForVectors
	disp('Input:');
	disp('------');
	vecA = [0.0, 0.3, -1.5, 0.0];
	disp('vecA = ');
	disp(vecA);
	vecB = [1.0, 0.2, 0.0, 1.0];
	disp('vecB = ');
	disp(vecB);
	disp('Test:');
	disp('-----');
	disp('vecA = vecA'':');
	vecA = vecA';
	disp('vecB = vecB'':');
	vecB = vecB';
	disp('B = SimilarityMatrixForVectors (vecA, vecB):');
	B = SimilarityMatrixForVectors (vecA, vecB);
	disp('B = ');
	disp(B);
	disp('err = norm(B*vecA-vecB):');
	err = norm(B*vecA-vecB);
	disp('err = ');
	disp(err);
	assert(norm(B*vecA-vecB)<10^-14 , 'The resulting matrix T does not satisfy T*vecA = vecB!');
end
