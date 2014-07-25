format short g;

disp('---BuTools: Discrete PH function test file---');

disp('Enable the verbose messages with the BuToolsVerbose flag')
global BuToolsVerbose;
BuToolsVerbose = true

disp('Enable input parameter checking with the BuToolsCheckInput flag')
global BuToolsCheckInput;
BuToolsCheckInput = true

disp('----------------------------------------------------------------------------');
help MomentsFromMG

disp('Input:');
disp('------');

a=[-0.6 0.3 1.3]
A=[0.25 0.2 -0.15; 0.3 0.1 0.25; 0 0.2 0.47]

disp('Test:');
disp('-----');

disp('MomentsFromMG(a,A,3):');
m=MomentsFromMG(a,A,3);
disp(m);

assert(all(m>0) && CheckMoments(m), 'MomentsFromMG returned invalid moments!');

disp('MomentsFromMG(a,A):');
m=MomentsFromMG(a,A);
disp(m);

assert(all(m>0) && CheckMoments(m), 'MomentsFromMG returned invalid moments!');

disp('----------------------------------------------------------------------------');
help MomentsFromDPH

disp('Input:');
disp('------');

a=[0.76 0 0.24]
A=[0.34 0.66 0; 0.79 0.05 0.07; 0.26 0.73 0.01]

disp('Test:');
disp('-----');

disp('MomentsFromDPH(a,A,5):');
m=MomentsFromDPH(a,A,5);
disp(m);

assert(all(m>0) && CheckMoments(m), 'MomentsFromDPH returned invalid moments!');

disp('----------------------------------------------------------------------------');
help PmfFromMG

disp('Input:');
disp('------');

a=[-0.6 0.3 1.3]
A=[0.25 0.2 -0.15; 0.3 0.1 0.25; 0 0.2 0.47]
x = (0:1:100)';

disp('Test:');
disp('-----');

disp('pmf=PmfFromMG(a,A,x):');
pmf = PmfFromMG(a,A,x);

assert(all(pmf)>=0, 'PmfFromMG returned negative pmf!');
assert(abs(pmf'*x - MomentsFromMG(a,A,1))<1e-12, 'The mean computed from the pmf does not match the theoretical result!');

disp('----------------------------------------------------------------------------');
help PmfFromDPH

disp('Input:');
disp('------');

a=[0.76 0 0.24]
A=[0.34 0.66 0; 0.79 0.05 0.07; 0.26 0.73 0.01]
x = (0:1:1000)';

disp('Test:');
disp('-----');

disp('pmf=PmfFromDPH(a,A,x):');
pmf = PmfFromDPH(a,A,x);

assert(all(pmf)>=0, 'PmfFromDPH returned negative pmf!');
assert(abs(pmf'*x - MomentsFromDPH(a,A,1))<1e-12, 'The mean computed from the pmf does not match the theoretical result!');

disp('----------------------------------------------------------------------------');
help CdfFromMG

disp('Input:');
disp('------');

a=[-0.6 0.3 1.3]
A=[0.25 0.2 -0.15; 0.3 0.1 0.25; 0 0.2 0.47]
x = (0:1:100)';

disp('Test:');
disp('-----');

disp('cdf=CdfFromMG(a,A,x):');
cdf = CdfFromMG(a,A,x);

assert(all(diff(cdf)>=0), 'The cdf is not increasing monotonously!');
assert(abs(sum(1-cdf) - MomentsFromMG(a,A,1))<1e-12, 'The mean computed from the cdf does not match the theoretical result!');

disp('----------------------------------------------------------------------------');
help CdfFromDPH

disp('Input:');
disp('------');

a=[0.76 0 0.24]
A=[0.34 0.66 0; 0.79 0.05 0.07; 0.26 0.73 0.01]
x = (0:1:1000)';

disp('Test:');
disp('-----');

disp('cdf=CdfFromDPH(a,A,x):');
cdf = CdfFromDPH(a,A,x);

assert(all(diff(cdf)>=0), 'The cdf is not increasing monotonously!');
assert(abs(sum(1-cdf) - MomentsFromDPH(a,A,1))<1e-12, 'The mean computed from the cdf does not match the theoretical result!');

disp('----------------------------------------------------------------------------');
help RandomDPH

disp('Test:');
disp('-----');
   
disp('[a,A]=RandomDPH(3,10,5):');
[a,A]=RandomDPH(3,10,5);
disp(a);
disp(A);

assert(CheckDPHRepresentation(a,A), 'RandomDPH failed to return a valid DPH representation!');
if max(abs(MomentsFromDPH(a,A,1)-10))>1e-14
    disp('RandomDPH failed to match the given mean value!');
end
clo = 1-sum(A,2);
assert (sum(a==0)+sum(sum(A==0))+sum(clo==0)==5, 'The number of zero entries does not match the function parameter!');

disp('----------------------------------------------------------------------------');
help CheckMGRepresentation

disp('Input:');
disp('------');

a=[-0.6 0.3 1.3]
A=[0.25 0.2 -0.15; 0.3 0.1 0.25; 0 0.2 0.47]

disp('Test:');
disp('-----');

disp('CheckMGRepresentation(a,A):');
flag=CheckMGRepresentation(a,A);
disp(flag);

assert(flag==1, 'CheckMGRepresentation failed to recognize a valid MG distribution!');

disp('Input:');
disp('------');

a=[-0.6 0.3 1.3]
A=[0.35 0.2 -0.25; 0.3 0.1 0.25; 0 0.2 0.47]

disp('Test:');
disp('-----');

disp('CheckMGRepresentation(a,A):');
flag=CheckMGRepresentation(a,A);
disp(flag);

assert(flag==0, 'CheckMGRepresentation failed to recognize wrong eigenvalues!');

disp('----------------------------------------------------------------------------');
help CheckDPHRepresentation

disp('Input:');
disp('------');

a=[0.48 0.08 0.26 0.18]
A=[0 0.08 0.08 0.8; 0.55 0 0.24 0.19; 0.06 0.03 0 0.001; 0.23 0.005 0.2 0.53]

disp('Test:');
disp('-----');

disp('CheckDPHRepresentation(a,A):');
flag=CheckDPHRepresentation(a,A);
disp(flag);

assert(flag==1, 'CheckDPHRepresentation failed to recognize a valid MG distribution!');

disp('Input:');
disp('------');

a=[0.48 0.08]
A=[0 0.08; 0.55 0.5]

disp('Test:');
disp('-----');

disp('CheckDPHRepresentation(a,A):');
flag=CheckDPHRepresentation(a,A);
disp(flag);

assert(flag==0, 'CheckDPHRepresentation failed to recognize wrong row sums!');

disp('----------------------------------------------------------------------------');
help MGFromMoments

disp('Input:');
disp('------');

moms = [4.08 20.41 130.45 1054.41 10463.73]

disp('Test:');
disp('-----');

disp('[b,B]=MGFromMoments(moms):');
[b,B]=MGFromMoments(moms)

disp('MomentsFromMG(b,B):');
m=MomentsFromMG(b,B);
disp(m);

assert(norm(moms-m)<1e-9,'The moments of the result of MGFromMoments do not match the input!');

disp('----------------------------------------------------------------------------');
help DPHFromMG

disp('Input:');
disp('------');

a=[-0.6 0.3 1.3]
A=[0.1 0.2 0; 0.3 0.1 0.25; -0.3 0.2 0.77]

disp('Test:');
disp('-----');

disp('CheckMGRepresentation(a,A):');
flag=CheckMGRepresentation(a,A);
disp(flag);

disp('CheckDPHRepresentation(a,A):');
flag=CheckDPHRepresentation(a,A);
disp(flag);

disp('[b,B]=DPHFromMG(a,A)');
[b,B]=DPHFromMG(a,A)

disp('Check the obtained DPH, CheckDPHRepresentation(b,B):');
flag=CheckDPHRepresentation(b,B);
disp(flag);

C=SimilarityMatrix(A,B);
err1 = max(max(abs(A*C-C*B)));
err2 = max(abs(a*C-b));

assert(flag && err1<1e-12 && err2<1e-12, 'Transformation to DPH failed!');


disp('----------------------------------------------------------------------------');
help CanonicalFromDPH2

disp('Input:');
disp('------');

a=[0 1]
A=[0.23 0.22; 0.41 0.48]

disp('Test:');
disp('-----');

disp('[b,B]=CanonicalFromDPH2(a,A):');
[b,B]=CanonicalFromDPH2(a,A)

disp('Eigenvalues of A:');
ev=EigSort(eig(A));
disp(ev);

disp('Check the obtained DPH, CheckDPHRepresentation(b,B):');
flag=CheckDPHRepresentation(b,B);
disp(flag);

C=SimilarityMatrix(A,B);
err1 = max(max(abs(A*C-C*B)));
err2 = max(abs(a*C-b));

assert(flag && err1<1e-12 && err2<1e-12, 'Transformation failed!');

disp('Input:');
disp('------');

a=[1 0]
A=[0 0.61; 0.56 0.44]

disp('Test:');
disp('-----');

disp('[b,B]=CanonicalFromDPH2(a,A):');
[b,B]=CanonicalFromDPH2(a,A)

disp('Eigenvalues of A:');
ev=EigSort(eig(A));
disp(ev);

disp('Check the obtained DPH, CheckDPHRepresentation(b,B):');
flag=CheckDPHRepresentation(b,B);
disp(flag);

C=SimilarityMatrix(A,B);
err1 = max(max(abs(A*C-C*B)));
err2 = max(abs(a*C-b));

assert(flag && err1<1e-12 && err2<1e-12, 'Transformation failed!');

disp('----------------------------------------------------------------------------');
help CanonicalFromDPH3

disp('Input:');
disp('------');

a=[0.46 0.22 0.32]
A=[0.67 0.01 0.12; 0.06 0.45 0.15; 0.18 0.43 0.32]

disp('Test:');
disp('-----');

disp('[b,B]=CanonicalFromDPH3(a,A):');
[b,B]=CanonicalFromDPH3(a,A)

disp('Eigenvalues of A:');
disp(EigSort(eig(A)));

disp('Check the obtained DPH, CheckDPHRepresentation(b,B):');
flag=CheckDPHRepresentation(b,B);
disp(flag);

C=SimilarityMatrix(A,B);
err1 = max(max(abs(A*C-C*B)));
err2 = max(abs(a*C-b));

assert(flag && err1<1e-12 && err2<1e-12, 'Transformation failed!');

disp('Input:');
disp('------');

a=[0.76 0.12 0.12]
A=[0.31 0 0; 0.98 0 0.02; 0.88 0.04 0.08]

disp('Test:');
disp('-----');

disp('[b,B]=CanonicalFromDPH3(a,A):');
[b,B]=CanonicalFromDPH3(a,A)

disp('Eigenvalues of A:');
disp(EigSort(eig(A)));

disp('Check the obtained DPH, CheckDPHRepresentation(b,B):');
flag=CheckDPHRepresentation(b,B);
disp(flag);

C=SimilarityMatrix(A,B);
err1 = max(max(abs(A*C-C*B)));
err2 = max(abs(a*C-b));

assert(flag && err1<1e-12 && err2<1e-12, 'Transformation failed!');

disp('Input:');
disp('------');

a=[0.67 0.07 0.26]
A=[0.31 0 0; 0.98 0 0.02; 0.88 0.04 0.08]

disp('Test:');
disp('-----');

disp('[b,B]=CanonicalFromDPH3(a,A):');
[b,B]=CanonicalFromDPH3(a,A)

disp('Eigenvalues of A:');
disp(EigSort(eig(A)));

disp('Check the obtained DPH, CheckDPHRepresentation(b,B):');
flag=CheckDPHRepresentation(b,B);
disp(flag);

C=SimilarityMatrix(A,B);
err1 = max(max(abs(A*C-C*B)));
err2 = max(abs(a*C-b));

assert(flag && err1<1e-12 && err2<1e-12, 'Transformation failed!');

disp('Input:');
disp('------');

a=[0.78 0.04 0.18]
A=[0.06 0.25 0.31; 0.45 0.18 0.33; 0.98 0 0.01]

disp('Test:');
disp('-----');

disp('[b,B]=CanonicalFromDPH3(a,A):');
[b,B]=CanonicalFromDPH3(a,A)

disp('Eigenvalues of A:');
disp(EigSort(eig(A)));

disp('Check the obtained DPH, CheckDPHRepresentation(b,B):');
flag=CheckDPHRepresentation(b,B);
disp(flag);

C=SimilarityMatrix(A,B);
err1 = max(max(abs(A*C-C*B)));
err2 = max(abs(a*C-b));

assert(flag && err1<1e-12 && err2<1e-12, 'Transformation failed!');

disp('----------------------------------------------------------------------------');
help AcyclicDPHFromMG

disp('Input:');
disp('------');

a=[0 0 1]
A=[0.22 0 0; 0.3 0.1 0.55; 0.26 0 0.73]

disp('Test:');
disp('-----');

disp('[b,B]=AcyclicDPHFromMG(a,A):');
[b,B]=AcyclicDPHFromMG(a,A)

disp('MomentsFromDPH(a,A),MomentsFromDPH(b,B):');
ma=MomentsFromDPH(a,A);
mb=MomentsFromDPH(b,B);
disp(ma);
disp(mb);

disp('Check the obtained DPH, CheckDPHRepresentation(b,B):');
flag=CheckDPHRepresentation(b,B);
disp(flag);

C=SimilarityMatrix(A,B);
err1 = max(max(abs(A*C-C*B)));
err2 = max(abs(a*C-b));

assert(flag && err1<1e-12 && err2<1e-12, 'Transformation failed!');

disp('----------------------------------------------------------------------------');
help DPH2From3Moments

disp('Input:');
disp('------');

a=[0.9 0.1];
A=[0.2 0.61; 0.58 0.41];
moms=MomentsFromDPH(a,A)

disp('Test:');
disp('-----');

disp('[b,B]=DPH2From3Moments(moms):');
[b,B]=DPH2From3Moments(moms)

disp('MomentsFromMG(b,B):');
m=MomentsFromMG(b,B);
disp(m);

assert(norm(moms-m)<1e-9,'The moments of the result of DPH2From3Moments do not match the input!');

disp('----------------------------------------------------------------------------');
help DPH3From5Moments

disp('Input:');
disp('------');

a=[0.7 0.1 0.2];
A=[0.2 0.51 0.1; 0.58 0.41 0; 0.1 0.4 0.3];
moms=MomentsFromDPH(a,A)

disp('Test:');
disp('-----');

disp('[b,B]=DPH3From5Moments(moms):');
[b,B]=DPH3From5Moments(moms)

disp('MomentsFromMG(b,B):');
m=MomentsFromMG(b,B);
disp(m);

assert(norm(moms-m)<1e-6,'The moments of the result of DPH3From5Moments do not match the input!');

disp('----------------------------------------------------------------------------');
help SamplesFromDPH

disp('Input:');
disp('------');

a=[0.76 0 0.24]
A=[0.34 0.66 0; 0.79 0.05 0.07; 0.26 0.73 0.01]

disp('Test:');
disp('-----');

disp('x=SamplesFromDPH(a,A,1000)');
x=SamplesFromDPH(a,A,1000);

disp('Moments from the samples:');
mt = MarginalMomentsFromTrace(x,3);
disp(mt);

disp('Moments from the DPH:');
mp = MomentsFromDPH(a,A,3);
disp(mp);
