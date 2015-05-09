format short g;

disp('---BuTools: Continous MAP function test file---');

disp('Enable the verbose messages with the BuToolsVerbose flag')
global BuToolsVerbose;
BuToolsVerbose = true

disp('Enable input parameter checking with the BuToolsCheckInput flag')
global BuToolsCheckInput;
BuToolsCheckInput = true

disp('----------------------------------------------------------------------------');
help MarginalDistributionFromMAP

disp('Input:');
disp('------');

D0=[-0.17 0 0 0.07; 0.01 -0.78 0.03 0.08; 0.22 0.17 -1.1 0.02; 0.04 0.12 0 -0.42]
D1=[0 0.06 0 0.04; 0.04 0.19 0.21 0.22; 0.22 0.13 0.15 0.19; 0.05 0 0.17 0.04]

disp('Test:');
disp('-----');

disp('[a,A]=MarginalDistributionFromMAP(D0,D1):');
[a,A]=MarginalDistributionFromMAP(D0,D1);
disp(a);
disp(A);

assert(length(a)==size(D0,2) && CheckPHRepresentation(a,A), 'MarginalDistributionFromMAP returned a wrong PH representation!');

disp('----------------------------------------------------------------------------');
help MarginalMomentsFromMAP

disp('Input:');
disp('------');

D0=[-0.17 0 0 0.07; 0.01 -0.78 0.03 0.08; 0.22 0.17 -1.1 0.02; 0.04 0.12 0 -0.42]
D1=[0 0.06 0 0.04; 0.04 0.19 0.21 0.22; 0.22 0.13 0.15 0.19; 0.05 0 0.17 0.04]

disp('Test:');
disp('-----');

disp('moms=MarginalMomentsFromMAP(D0,D1):');
moms=MarginalMomentsFromMAP(D0,D1);
disp(moms);

assert(length(moms)==2*size(D0,2)-1 && CheckMoments(moms), 'MarginalMomentsFromMAP returned wrong moments!');

disp('----------------------------------------------------------------------------');
help MarginalDistributionFromRAP

disp('Input:');
disp('------');

H0=[-2 0 0; 0 -3 1; 0 -1 -2];
H1=[1.8 0.2 0; 0.2 1.8 0; 0.2 1.8 1];

disp('Test:');
disp('-----');

disp('[a,A]=MarginalDistributionFromRAP(H0,H1):');
[a,A]=MarginalDistributionFromRAP(H0,H1);
disp(a);
disp(A);

assert(length(a)==size(H0,2) && CheckMERepresentation(a,A), 'MarginalDistributionFromRAP returned a wrong ME representation!');

disp('----------------------------------------------------------------------------');
help MarginalMomentsFromRAP

disp('Input:');
disp('------');

H0=[-2 0 0; 0 -3 1; 0 -1 -2];
H1=[1.8 0.2 0; 0.2 1.8 0; 0.2 1.8 1];

disp('Test:');
disp('-----');

disp('moms=MarginalMomentsFromRAP(H0,H1):');
moms=MarginalMomentsFromRAP(H0,H1);
disp(moms);

assert(length(moms)==2*size(H0,2)-1 && CheckMoments(moms), 'MarginalMomentsFromRAP returned wrong moments!');

disp('----------------------------------------------------------------------------');
help MarginalDistributionFromMMAP

D0=[-1.78 0.29; 0.07 -0.92]
D1=[0.15 0.49; 0.23 0.36]
D2=[0.11 0.2; 0.01 0]
D3=[0.14 0.4; 0.11 0.14]

disp('Test:');
disp('-----');

disp('[a,A]=MarginalDistributionFromMMAP({D0,D1,D2,D3}):');
[a,A]=MarginalDistributionFromMMAP({D0,D1,D2,D3});
disp(a);
disp(A);

assert(length(a)==size(D0,2) && CheckPHRepresentation(a,A), 'MarginalDistributionFromMMAP returned a wrong PH representation!');

disp('----------------------------------------------------------------------------');
help MarginalMomentsFromMMAP

D0=[-1.78 0.29; 0.07 -0.92]
D1=[0.15 0.49; 0.23 0.36]
D2=[0.11 0.2; 0.01 0]
D3=[0.14 0.4; 0.11 0.14]

disp('Test:');
disp('-----');

disp('moms=MarginalMomentsFromMMAP({D0,D1,D2,D3}):');
moms=MarginalMomentsFromMMAP({D0,D1,D2,D3});
disp(moms);

assert(length(moms)==2*size(D0,2)-1 && CheckMoments(moms), 'MarginalMomentsFromMMAP returned wrong moments!');

disp('----------------------------------------------------------------------------');
help MarginalDistributionFromMRAP

disp('Input:');
disp('------');

x=0.18;
H0=[-5 0.1+x 0.9 1; 1 -8 0.9 0.1; 0.9 0.1 -4 1; 1 2 3 -9]
H1=[0.1-x 0.7 0.1 0.1; 0.1 1 1.8 0.1; 0.1 0.1 0.1 0.7; 0.7 0.1 0.1 0.1]
H2=[0.1 0.1 0.1 1.7; 1.8 0.1 1 0.1; 0.1 0.1 0.7 0.1; 0.1 1 0.1 0.8]

disp('Test:');
disp('-----');

disp('[a,A]=MarginalDistributionFromMRAP({H0,H1,H2}):');
[a,A]=MarginalDistributionFromMRAP({H0,H1,H2});
disp(a);
disp(A);

assert(length(a)==size(H0,2) && CheckMERepresentation(a,A), 'MarginalDistributionFromMRAP returned a wrong ME representation!');

disp('----------------------------------------------------------------------------');
help MarginalMomentsFromMRAP

disp('Input:');
disp('------');

x=0.18;
H0=[-5 0.1+x 0.9 1; 1 -8 0.9 0.1; 0.9 0.1 -4 1; 1 2 3 -9]
H1=[0.1-x 0.7 0.1 0.1; 0.1 1 1.8 0.1; 0.1 0.1 0.1 0.7; 0.7 0.1 0.1 0.1]
H2=[0.1 0.1 0.1 1.7; 1.8 0.1 1 0.1; 0.1 0.1 0.7 0.1; 0.1 1 0.1 0.8]

disp('Test:');
disp('-----');

disp('moms=MarginalMomentsFromMRAP({H0,H1,H2}):');
moms=MarginalMomentsFromMRAP({H0,H1,H2});
disp(moms);

assert(length(moms)==2*size(H0,2)-1 && CheckMoments(moms), 'MarginalMomentsFromMRAP returned wrong moments!');

disp('----------------------------------------------------------------------------');
help LagCorrelationsFromMAP

disp('Input:');
disp('------');

D0=[-5 0 1 1; 1 -8 1 0; 1 0 -4 1; 1 2 3 -9]
D1=[0 1 0 2; 2 1 3 0; 0 0 1 1; 1 1 0 1]

disp('Test:');
disp('-----');

disp('LagCorrelationsFromMAP(D0,D1,3):');
corr = LagCorrelationsFromMAP(D0,D1,3);
disp(corr);

assert(length(corr)==3 && all(corr<1) && all(corr>-1), 'LagCorrelationsFromMAP returned wrong autocorrelation coefficients!');

disp('----------------------------------------------------------------------------');
help LagCorrelationsFromRAP

disp('Input:');
disp('------');

H0=[-2 0 0; 0 -3 1; 0 -1 -2]
H1=[1.8 0.2 0; 0.2 1.8 0; 0.2 1.8 1]

disp('Test:');
disp('-----');

disp('LagCorrelationsFromRAP(H0,H1,3):');
corr = LagCorrelationsFromRAP(H0,H1,3);
disp(corr);

assert(length(corr)==3 && all(corr<1) && all(corr>-1), 'LagCorrelationsFromRAP returned wrong autocorrelation coefficients!');

disp('----------------------------------------------------------------------------');
help LagkJointMomentsFromMAP

disp('Input:');
disp('------');

D0=[-5 0 1 1; 1 -8 1 0; 1 0 -4 1; 1 2 3 -9]
D1=[0 1 0 2; 2 1 3 0; 0 0 1 1; 1 1 0 1]

disp('Test:');
disp('-----');

disp('LagkJointMomentsFromMAP(D0,D1,2,1):');
Nm=LagkJointMomentsFromMAP(D0,D1,4,1);
disp(Nm);

moms=MarginalMomentsFromMAP(D0,D1,4);

disp('Correlation from joint moments:');
cjm=zeros(3,1);
for i=1:3
    Nx=LagkJointMomentsFromMAP(D0,D1,1,i);
    cjm(i) = (Nx(2,2)-moms(1)^2) / (moms(2)-moms(1)^2);
end
disp(cjm)
corr = LagCorrelationsFromMAP(D0,D1,3);

assert(all(all(Nm>0)) && norm(moms-Nm(1,2:end))<1e-14 && norm(moms'-Nm(2:end,1))<1e-14 && norm(corr-cjm)<1e-14, 'Joint moment matrix is invalid!');

disp('----------------------------------------------------------------------------');
help LagkJointMomentsFromRAP

disp('Input:');
disp('------');

H0=[-2 0 0; 0 -3 1; 0 -1 -2];
H1=[1.8 0.2 0; 0.2 1.8 0; 0.2 1.8 1];

disp('Test:');
disp('-----');

disp('LagkJointMomentsFromRAP(H0,H1,2,1):');
Nm=LagkJointMomentsFromRAP(H0,H1,4,1);
disp(Nm);

moms=MarginalMomentsFromRAP(H0,H1,4);

disp('Correlation from joint moments:');
cjm=zeros(3,1);
for i=1:3
    Nx=LagkJointMomentsFromRAP(H0,H1,1,i);
    cjm(i) = (Nx(2,2)-moms(1)^2) / (moms(2)-moms(1)^2);
end
disp(cjm)
corr = LagCorrelationsFromRAP(H0,H1,3);

assert(all(all(Nm>0)) && norm(moms-Nm(1,2:end))<1e-14 && norm(moms'-Nm(2:end,1))<1e-14 && norm(corr-cjm)<1e-14, 'Joint moment matrix is invalid!');

disp('----------------------------------------------------------------------------');
help LagkJointMomentsFromMMAP

disp('Input:');
disp('------');

D0=[-1.78 0.29; 0.07 -0.92]
D1=[0.15 0.49; 0.23 0.36]
D2=[0.11 0.2; 0.01 0]
D3=[0.14 0.4; 0.11 0.14]

disp('Test:');
disp('-----');

disp('LagkJointMomentsFromMMAP({D0,D1,D2,D3},3,1)');
Nm=LagkJointMomentsFromMMAP({D0,D1,D2,D3},3,1);

disp('Moments of arrival type 1, Nm{1}:');
disp(Nm{1});
disp('Moments of arrival type 2, Nm{2}:');
disp(Nm{2});
disp('Moments of arrival type 2, Nm{3}:');
disp(Nm{3});

assert(length(Nm)==3 && norm(SumMatrixList(Nm)-LagkJointMomentsFromMAP(D0,D1+D2+D3,3,1))<1e-14,'Joint moment matrix is invalid!');

disp('----------------------------------------------------------------------------');
help LagkJointMomentsFromMRAP

disp('Input:');
disp('------');

x=0.18;
H0=[-5 0.1+x 0.9 1; 1 -8 0.9 0.1; 0.9 0.1 -4 1; 1 2 3 -9]
H1=[0.1-x 0.7 0.1 0.1; 0.1 1 1.8 0.1; 0.1 0.1 0.1 0.7; 0.7 0.1 0.1 0.1]
H2=[0.1 0.1 0.1 1.7; 1.8 0.1 1 0.1; 0.1 0.1 0.7 0.1; 0.1 1 0.1 0.8]
disp('Test:');
disp('-----');

disp('LagkJointMomentsFromMRAP({H0,H1,H2},3,2)');
Nm=LagkJointMomentsFromMRAP({H0,H1,H2},3,2);

disp('Moments of arrival type 1, Nm{1}:');
disp(Nm{1});
disp('Moments of arrival type 2, Nm{2}:');
disp(Nm{2});

assert(length(Nm)==2 && norm(SumMatrixList(Nm)-LagkJointMomentsFromRAP(H0,H1+H2,3,2))<1e-14,'Joint moment matrix is invalid!');

disp('----------------------------------------------------------------------------');
help RandomMAP

disp('Test:');
disp('-----');

disp('[D0,D1]=RandomMAP(4,1.62,10):');
[D0,D1]=RandomMAP(4,1.62,10);
disp(D0);
disp(D1);

disp('Check the mean of the obtained MAP:');
m = MarginalMomentsFromMAP(D0,D1,1);
disp(m);

assert(CheckMAPRepresentation(D0,D1), 'RandomMAP failed to return a valid MAP representation!');
assert(max(abs(m-1.62))<1e-14, 'RandomMAP failed to match the given mean value!');

disp('----------------------------------------------------------------------------');
help RandomMMAP

disp('Test:');
disp('-----');

disp('D=RandomMMAP(4,3,1.62,10):');
D=RandomMMAP(4,3,1.62,10);
disp(D{1});
disp(D{2});
disp(D{3});
disp(D{4});

disp('Check the mean of the obtained MMAP:');
m = MarginalMomentsFromMMAP(D,1);
disp(m);

assert(CheckMMAPRepresentation(D), 'RandomMMAP failed to return a valid MMAP representation!');
assert(max(abs(m-1.62))<1e-14, 'RandomMMAP failed to match the given mean value!');

disp('----------------------------------------------------------------------------');
help CheckMAPRepresentation

disp('Input:');
disp('------');

D0=[-1 0 1; 0 -2 0; 1 0 -3]
D1=[-1 0 1 0; 0 -2 0 1; 1 0 -3 0; 1 2 2 1]

disp('Test:');
disp('-----');

disp('CheckMAPRepresentation(D0,D1):');
flag=CheckMAPRepresentation(D0,D1);
disp(flag);

assert(flag==0, 'CheckMAPRepresentation failed to detect the incompatible shapes of D0 and D1!');

disp('Input:');
disp('------');

D0=[-1 0 1; 0 -2 0; 1 0 -3]
D1=[1 0 1; 0 2 0 ; 1 0 3]

disp('Test:');
disp('-----');

disp('CheckMAPRepresentation(D0,D1):');
flag=CheckMAPRepresentation(D0,D1);
disp(flag);

assert(flag==0, 'CheckMAPRepresentation failed to detect invalid rowsums!');

disp('Input:');
disp('------');

D0=[-3 0 1; 0 -2 0; 1 0 -5]
D1=[1 0 1; 0 2 0 ; 1 0 3]

disp('Test:');
disp('-----');

disp('CheckMAPRepresentation(D0,D1):');
flag=CheckMAPRepresentation(D0,D1);
disp(flag);

assert(flag==1, 'CheckMAPRepresentation failed to recognize a valid MAP representation!');

disp('----------------------------------------------------------------------------');
help CheckRAPRepresentation

disp('Input:');
disp('------');

H0=[-1 0 1; 0 -2 0; 1 0 -3; 1 2 2]
H1=[-1 0 1; 0 -2 0; 1 0 -3; 1 2 2]

disp('Test:');
disp('-----');

disp('CheckRAPRepresentation(H0,H1):')
flag=CheckRAPRepresentation(H0,H1);
disp(flag);

assert(flag==0, 'CheckRAPRepresentation failed to detect the incompatible shapes of D0 and D1!');

disp('Input:');
disp('------');

H0=[-1 0 2; 0 2 0; 1 0 -3]
H1=[-1 0 1; 0 -2 0; 1 0 -3]

disp('Test:');
disp('-----');

disp('CheckRAPRepresentation(H0,H1):')
flag=CheckRAPRepresentation(H0,H1);
disp(flag);

assert(flag==0, 'CheckRAPRepresentation failed to detect invalid rowsums!');

disp('Input:');
disp('------');

H0=[-1 0 0; 0 -2 2; 0 3 -3]
H1=[0 0 1; 0 -1 1; 1 0 -1]

disp('Test:');
disp('-----');

disp('CheckRAPRepresentation(H0,H1):')
flag=CheckRAPRepresentation(H0,H1);
disp(flag);

assert(flag==0, 'CheckRAPRepresentation failed to detect invalid eigenvalues!');

disp('Input:');
disp('------');

H0=[-2 0 0; 0 -1 1; 0 -1 -1]
H1=[1 0 1; 0 1 -1; 1 0 1]

disp('Test:');
disp('-----');

disp('CheckRAPRepresentation(H0,H1):')
flag=CheckRAPRepresentation(H0,H1);
disp(flag);

assert(flag==0, 'CheckRAPRepresentation failed to detect invalid eigenvalues!');

disp('Input:');
disp('------');

H0=[-1 0 0; 0 -2 1; 0 -1 -2]
H1=[1 0 0; 0 1 0; 1 1 1]

disp('Test:');
disp('-----');

disp('CheckRAPRepresentation(H0,H1):')
flag=CheckRAPRepresentation(H0,H1);
disp(flag);

assert(flag==1, 'CheckRAPRepresentation failed to recognize a valid RAP representation!');

disp('----------------------------------------------------------------------------');
help CheckMMAPRepresentation

disp('Input:');
disp('------');

D0=[-1.05 0.03 0.07; 0.19 -1.63 0.06; 0 0.2 -1.03]
D1=[0.16 0.11 0; 0.1 0.16 0; 0.27 0 0.19]
D2=[0.01 0.09 0.13; 0.26 0.21 0.05; 0 0.16 0.07]
D3=[0.19 0.06 0.2; 0.17 0.16 0.27; 0 0 0.14]

disp('Test:');
disp('-----');

disp('CheckMMAPRepresentation({D0,D1,D2,D3}):');
flag=CheckMMAPRepresentation({D0,D1,D2,D3});
disp(flag);

assert(flag==1, 'CheckMMAPRepresentation failed to recognize a valid MMAP representation!');

disp('----------------------------------------------------------------------------');
help CheckMRAPRepresentation

disp('Input:');
disp('------');

x=0.18;
H0=[-5 0.1+x 0.9 1; 1 -8 0.9 0.1; 0.9 0.1 -4 1; 1 2 3 -9]
H1=[0.1-x 0.7 0.1 0.1; 0.1 1 1.8 0.1; 0.1 0.1 0.1 0.7; 0.7 0.1 0.1 0.1]
H2=[0.1 0.1 0.1 1.7; 1.8 0.1 1 0.1; 0.1 0.1 0.7 0.1; 0.1 1 0.1 0.8]

disp('Test:');
disp('-----');

disp('CheckMRAPRepresentation({H0,H1,H2}):');
flag=CheckMRAPRepresentation({H0,H1,H2});
disp(flag);

assert(flag==1, 'CheckMRAPRepresentation failed to recognize a valid MRAP representation!');

disp('----------------------------------------------------------------------------');
help RAPFromMoments

disp('Input:');
disp('------');

G0=[-6.2 2 0; 2 -9 1; 1 0 -3];
G1=[2.2 -2 4; 2 2 2; 1 0 1];
moms=MarginalMomentsFromRAP(G0,G1,5)
Nm=LagkJointMomentsFromRAP(G0,G1,2,1)

disp('Test:');
disp('-----');

disp('[H0,H1]=RAPFromMoments(moms,Nm):');
[H0,H1]=RAPFromMoments(moms,Nm);
disp(H0);
disp(H1);

rmoms=MarginalMomentsFromRAP(H0,H1,5)
rNm=LagkJointMomentsFromRAP(H0,H1,2,1)

assert(norm(moms-rmoms)<1e-12 && norm(Nm-rNm)<1e-12, 'The moments and joint moments returned by RAPFromMoments are not the same as given!');

disp('Input:');
disp('------');

G0=[-5 0 1 1; 1 -8 1 0; 1 0 -4 1; 1 2 3 -9];
G1=[0 1 0 2; 2 1 3 0; 0 0 1 1; 1 1 0 1];
moms=MarginalMomentsFromRAP(G0,G1,7)
Nm=LagkJointMomentsFromRAP(G0,G1,3,1)

disp('Test:');
disp('-----');

disp('[H0,H1]=RAPFromMoments(moms,Nm):');
[H0,H1]=RAPFromMoments(moms,Nm);
disp(H0);
disp(H1);

global BuToolsCheckPrecision;
BuToolsCheckPrecision = 1e-8;
rmoms=MarginalMomentsFromRAP(H0,H1,7)
rNm=LagkJointMomentsFromRAP(H0,H1,3,1)

assert(CheckRAPRepresentation(H0,H1,1e-8), 'RAPFromMoments returned an invalid RAP representation!');
assert(norm(moms-rmoms)<1e-8 && norm(Nm-rNm)<1e-8, 'The moments and joint moments returned by RAPFromMoments are not the same as given!');

disp('----------------------------------------------------------------------------');
help MRAPFromMoments

disp('Input:');
disp('------');

G0=[-1.05 0.03 0.07; 0.19 -1.63 0.06; 0 0.2 -1.03]
G1=[0.16 0.11 0; 0.1 0.16 0; 0.27 0 0.19]
G2=[0.01 0.09 0.13; 0.26 0.21 0.05; 0 0.16 0.07]
G3=[0.19 0.06 0.2; 0.17 0.16 0.27; 0 0 0.14]
G={G0,G1,G2,G3};
moms=MarginalMomentsFromMRAP(G,5)
Nm=LagkJointMomentsFromMRAP(G,2,1);
Nm1=Nm{1}
Nm2=Nm{2}
Nm3=Nm{3}

disp('Test:');
disp('-----');

disp('H=MRAPFromMoments(moms,Nm):');
H=MRAPFromMoments(moms,Nm);
disp(H{1});
disp(H{2});
disp(H{3});
disp(H{4});

BuToolsCheckPrecision = 1e-11;
rmoms=MarginalMomentsFromMRAP(H,5)
rNm=LagkJointMomentsFromMRAP(H,2,1);
rNm1=rNm{1}
rNm2=rNm{2}
rNm3=rNm{3}

assert(norm(moms-rmoms)<1e-10 && norm(Nm1-rNm1)<1e-10 && norm(Nm2-rNm2)<1e-10 && norm(Nm3-rNm3)<1e-10, 'The moments and joint moments returned by MRAPFromMoments are not the same as given!');

disp('----------------------------------------------------------------------------');
help RAPFromMomentsAndCorrelations

disp('Input:');
disp('------');

H0=[-6.2 2 0; 2 -9 1; 1 0 -3];
H1=[2.2 0 2; 0 4 2; 0 1 1];
mom=MarginalMomentsFromRAP(H0,H1)
corr=LagCorrelationsFromRAP(H0,H1,3)

disp('Test:');
disp('-----');

disp('[G0,G1]=RAPFromMomentsAndCorrelations(mom,corr)');
[G0,G1]=RAPFromMomentsAndCorrelations(mom,corr);
disp(G0);
disp(G1);

rmom=MarginalMomentsFromRAP(G0,G1);
rcorr=LagCorrelationsFromRAP(G0,G1,3);

assert(CheckRAPRepresentation(G0,G1), 'RAPFromMomentsAndCorrelations returned an invalid RAP representation!');
assert(norm(rmom-mom)+norm(rcorr-corr)<1e-12, 'The result of RAPFromMomentsAndCorrelations has different moments or correlations than given!');

disp('----------------------------------------------------------------------------');
help MAP2FromMoments

disp('Input:');
disp('------');

D0=[-14 1; 1 -25];
D1=[6 7; 3 21];
moms=MarginalMomentsFromMAP(D0, D1, 3)
corr=LagCorrelationsFromMAP(D0, D1, 1)

disp('Test:');
disp('-----');

disp('[D0,D1]=MAP2FromMoments(moms,corr):');
[D0,D1]=MAP2FromMoments(moms,corr);
disp(D0);
disp(D1);

rmoms=MarginalMomentsFromMAP(D0, D1, 3)
rcorr=LagCorrelationsFromMAP(D0, D1, 1)

assert(CheckMAPRepresentation(D0,D1), 'MAP2FromMoments returned an invalid MAP representation!');
assert(norm(moms-rmoms)<1e-12 && norm(corr-rcorr)<1e-12, 'The moments and the correlation returned by MAP2FromMoments are not the same as given!');

disp('----------------------------------------------------------------------------');
help MAP2CorrelationBounds

disp('Input:');
disp('------');

D0=[-14 1; 1 -25];
D1=[6 7; 3 21];
moms=MarginalMomentsFromMAP(D0,D1,3)

disp('Test:');
disp('-----');

disp('[lb,ub]=MAP2CorrelationBounds(moms):');
[lb,ub]=MAP2CorrelationBounds(moms);
disp(lb)
disp(ub)

assert(lb<=0 && lb>=-1 && ub>=0 && ub<=1, 'Correlation bounds given by MAP2CorrelationBounds are not correct');

disp('----------------------------------------------------------------------------');
help MAPFromFewMomentsAndCorrelations

disp('Input:');
disp('------');

moms = [1.1, 6.05]
corr1 = -0.17

disp('Test:');
disp('-----');

disp('[D0,D1]=MAPFromFewMomentsAndCorrelations(moms, corr1):');
[D0,D1]=MAPFromFewMomentsAndCorrelations(moms, corr1);
disp(D0);
disp(D1);

rmoms = MarginalMomentsFromMAP(D0,D1,2);
rcorr1 = LagCorrelationsFromMAP(D0,D1,1);

assert(CheckMAPRepresentation(D0,D1), 'MAPFromFewMomentsAndCorrelations returned with a non-Markovian representation!');
assert(norm(rmoms-moms)<1e-12 && norm(rcorr1-corr1)<1e-12, 'MAPFromFewMomentsAndCorrelations failed to match the marginal moments or the lag-1 autocorrelation!');

disp('Input:');
disp('------');

moms = [1.2, 4.32, 20]
corr1 = -0.4

disp('Test:');
disp('-----');

disp('[D0,D1]=MAPFromFewMomentsAndCorrelations(moms, corr1):');
[D0,D1]=MAPFromFewMomentsAndCorrelations(moms, corr1);
disp(D0);
disp(D1);

rmoms = MarginalMomentsFromMAP(D0,D1,3);
rcorr1 = LagCorrelationsFromMAP(D0,D1,1);

assert(CheckMAPRepresentation(D0,D1,1e-13), 'MAPFromFewMomentsAndCorrelations returned with a non-Markovian representation!');
assert(norm(rmoms-moms)<1e-12 && norm(rcorr1-corr1)<1e-12, 'MAPFromFewMomentsAndCorrelations failed to match the marginal moments or the lag-1 autocorrelation!');

disp('Input:');
disp('------');

moms = [1.2, 4.32, 20]
corr1 = 0.4

disp('Test:');
disp('-----');

disp('[D0,D1]=MAPFromFewMomentsAndCorrelations(moms, corr1):');
[D0,D1]=MAPFromFewMomentsAndCorrelations(moms, corr1);
disp(D0);
disp(D1);

rmoms = MarginalMomentsFromMAP(D0,D1,3);
rcorr1 = LagCorrelationsFromMAP(D0,D1,1);

assert(CheckMAPRepresentation(D0,D1,1e-13), 'MAPFromFewMomentsAndCorrelations returned with a non-Markovian representation!');
assert(norm(rmoms-moms)<1e-12 && norm(rcorr1-corr1)<1e-12, 'MAPFromFewMomentsAndCorrelations failed to match the marginal moments or the lag-1 autocorrelation!');

disp('----------------------------------------------------------------------------');
help CanonicalFromMAP2

disp('Input:');
disp('------');

D0=[-14 1; 1 -25]
D1=[6 7; 3 21]

disp('Test:');
disp('-----');

disp('[H0,H1]=CanonicalFromMAP2(D0,D1):');
[H0,H1]=CanonicalFromMAP2(D0,D1);
disp(H0);
disp(H1);

C=SimilarityMatrix(H0,D0);
err = norm(H0*C-C*D0) + norm(H1*C-C*D1);

assert(CheckMAPRepresentation(H0,H1), 'The result of CanonicalFromMAP2 is not a valid MAP representation!');
assert(err<1e-12, 'The MAP returned by CanonicalFromMAP2 is not similar to the input!');

disp('----------------------------------------------------------------------------');
help MAPFromRAP

disp('Input:');
disp('------');

D0=[-2 2; 2 -9]
D1=[-2 2; 3 4]

disp('Test:');
disp('-----');

disp('[H0,H1]=MAPFromRAP(D0,D1):');
[H0,H1]=MAPFromRAP(D0,D1);
disp(H0);
disp(H1);

err = norm(LagkJointMomentsFromRAP(D0,D1,3,1)-LagkJointMomentsFromRAP(H0,H1,3,1));
assert(err<1e-12, 'The RAP returned by MAPFromRAP is not similar to the input!');

disp('Input:');
disp('------');

D0=[-2.4 2; 2 -9]
D1=[-1.6 2; 3 4]

disp('Test:');
disp('-----');

disp('[H0,H1]=MAPFromRAP(D0,D1):');
[H0,H1]=MAPFromRAP(D0,D1);
disp(H0);
disp(H1);

err = norm(LagkJointMomentsFromRAP(D0,D1,3,1)-LagkJointMomentsFromRAP(H0,H1,3,1));
assert(err<1e-12, 'The MAP returned by MAPFromRAP is not similar to the input!');
assert(CheckMAPRepresentation(H0,H1), 'The result of MAPFromRAP is not a MAP, as it should be!');

disp('----------------------------------------------------------------------------');
help MMAPFromMRAP

disp('Input:');
disp('------');

x=0.18;
H0=[-5 0.1+x 0.9 1; 1 -8 0.9 0.1; 0.9 0.1 -4 1; 1 2 3 -9]
H1=[0.1-x 0.7 0.1 0.1; 0.1 1 1.8 0.1; 0.1 0.1 0.1 0.7; 0.7 0.1 0.1 0.1]
H2=[0.1 0.1 0.1 1.7; 1.8 0.1 1 0.1; 0.1 0.1 0.7 0.1; 0.1 1 0.1 0.8]
H={H0,H1,H2};

moms=MarginalMomentsFromMRAP(H);
jmom=LagkJointMomentsFromMRAP(H,3,1);

disp('Test:');
disp('-----');

disp('G=MMAPFromMRAP({H0,H1,H2}):');
G=MMAPFromMRAP(H);
disp(G{1})
disp(G{2})
disp(G{3})

rmoms=MarginalMomentsFromMMAP(G);
rjmom=LagkJointMomentsFromMMAP(G,3,1);

err = norm(rjmom{1}-jmom{1}) + norm(rjmom{2}-jmom{2});
assert(err<1e-12, 'The MMAP returned by MMAPFromMRAP is not similar to the input!');
assert(CheckMMAPRepresentation(G), 'The result of MMAPFromMRAP is not a MMAP, as it should be!');

disp('----------------------------------------------------------------------------');
help MinimalRepFromRAP

disp('Input:');
disp('------');

D0=[-5 1 0; 3 -3 0; 1 1 -5]
D1=[0 0 4; 0 0 0; 1 1 1]

disp('Test:');
disp('-----');

disp('[H0,H1]=MinimalRepFromRAP(D0,D1,''cont'')');
[H0,H1]=MinimalRepFromRAP(D0,D1,'cont');
disp(H0);
disp(H1);

C = SimilarityMatrix(H0,D0);
err = norm(H0*C-C*D0) + norm(H1*C-C*D1);

assert(CheckRAPRepresentation(H0,H1), 'MinimalRepFromRAP did not return a valid RAP representation!');
assert(size(H0,1)==3 && err<1e-12, 'MinimalRepFromRAP returned a RAP which is non-similar to the input or has an unexpected size!');

disp('[H0,H1]=MinimalRepFromRAP(D0,D1,''obs'')');
[H0,H1]=MinimalRepFromRAP(D0,D1,'obs');
disp(H0);
disp(H1);

C = SimilarityMatrix(H0,D0);
err = norm(H0*C-C*D0) + norm(H1*C-C*D1);

assert(CheckRAPRepresentation(H0,H1), 'MinimalRepFromRAP did not return a valid RAP representation!');
assert(size(H0,1)==2 && err<1e-12, 'MinimalRepFromRAP returned a RAP which is non-similar to the input or has an unexpected size!');

disp('[H0,H1]=MinimalRepFromRAP(D0,D1,''obscont'')');
[H0,H1]=MinimalRepFromRAP(D0,D1,'obscont');
disp(H0);
disp(H1);

C = SimilarityMatrix(H0,D0);
err = norm(H0*C-C*D0) + norm(H1*C-C*D1);

assert(CheckRAPRepresentation(H0,H1), 'MinimalRepFromRAP did not return a valid RAP representation!');
assert(size(H0,1)==2 && err<1e-12, 'MinimalRepFromRAP returned a RAP which is non-similar to the input or has an unexpected size!');

disp('----------------------------------------------------------------------------');
help MinimalRepFromMRAP

disp('Input:');
disp('------');

D0=[-5 1 0; 3 -3 0; 1 1 -5]
D1=0.2*[0 0 4; 0 0 0; 1 1 1]
D2=0.8*[0 0 4; 0 0 0; 1 1 1]
D={D0,D1,D2};

disp('Test:');
disp('-----');

disp('H=MinimalRepFromMRAP(D,''cont'')');
H=MinimalRepFromMRAP(D,'cont');
disp(H{1});
disp(H{2});
disp(H{3});

C = SimilarityMatrix(H{1},D{1});
err = norm(H{1}*C-C*D{1}) + norm(H{2}*C-C*D{2}) + norm(H{3}*C-C*D{3});

assert(CheckMRAPRepresentation(H), 'MinimalRepFromMRAP did not return a valid MRAP representation!');
assert(size(H{1},1)==3 && err<1e-12, 'MinimalRepFromMRAP returned a MRAP which is non-similar to the input or has an unexpected size!');

disp('H=MinimalRepFromMRAP(D,''obs'')');
H=MinimalRepFromMRAP(D,'obs');
disp(H{1});
disp(H{2});
disp(H{3});

C = SimilarityMatrix(H{1},D{1});
err = norm(H{1}*C-C*D{1}) + norm(H{2}*C-C*D{2}) + norm(H{3}*C-C*D{3});

assert(CheckMRAPRepresentation(H), 'MinimalRepFromMRAP did not return a valid MRAP representation!');
assert(size(H{1},1)==2 && err<1e-12, 'MinimalRepFromMRAP returned a MRAP which is non-similar to the input or has an unexpected size!');

disp('H=MinimalRepFromMRAP(D,''obscont'')');
H=MinimalRepFromMRAP(D,'obscont');
disp(H{1});
disp(H{2});
disp(H{3});

C = SimilarityMatrix(H{1},D{1});
err = norm(H{1}*C-C*D{1}) + norm(H{2}*C-C*D{2}) + norm(H{3}*C-C*D{3});

assert(CheckMRAPRepresentation(H), 'MinimalRepFromMRAP did not return a valid MRAP representation!');
assert(size(H{1},1)==2 && err<1e-12, 'MinimalRepFromMRAP returned a MRAP which is non-similar to the input or has an unexpected size!');

disp('----------------------------------------------------------------------------');
help SamplesFromMAP

disp('Input:');
disp('------');

D0=[-0.17 0 0 0.07; 0.01 -0.78 0.03 0.08; 0.22 0.17 -1.1 0.02; 0.04 0.12 0 -0.42]
D1=[0 0.06 0 0.04; 0.04 0.19 0.21 0.22; 0.22 0.13 0.15 0.19; 0.05 0 0.17 0.04]

disp('Test:');
disp('-----');

disp('x=SamplesFromMAP(D0,D1,10000)');
x=SamplesFromMAP(D0,D1,10000);

disp('Moments from the samples:');
mt = MarginalMomentsFromTrace(x,3);
disp(mt);

disp('Moments from the MAP:');
mm = MarginalMomentsFromMAP(D0,D1,3);
disp(mm);

disp('----------------------------------------------------------------------------');
help SamplesFromMMAP

disp('Input:');
disp('------');

D0=[-1.78 0.29; 0.07 -0.92]
D1=[0.15 0.49; 0.23 0.36]
D2=[0.11 0.2; 0.01 0]
D3=[0.14 0.4; 0.11 0.14]
D = {D0,D1,D2,D3};

disp('Test:');
disp('-----');

disp('x=SamplesFromMMAP(D,10000)');
x=SamplesFromMMAP(D,10000);

disp('Moments from the samples:');
mt = MarginalMomentsFromTrace(x(:,1),3);
disp(mt);

disp('Moments from the MMAP:');
mm = MarginalMomentsFromMMAP(D,3);
disp(mm);
