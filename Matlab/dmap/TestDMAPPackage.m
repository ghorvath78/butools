format short g;

disp('---BuTools: Discrete MAP function test file---');

disp('Enable the verbose messages with the BuToolsVerbose flag')
global BuToolsVerbose;
BuToolsVerbose = true

disp('Enable input parameter checking with the BuToolsCheckInput flag')
global BuToolsCheckInput;
BuToolsCheckInput = true

disp('----------------------------------------------------------------------------');
help MarginalDistributionFromDMAP

disp('Input:');
disp('------');

D0=[0 0.02 0 0; 0 0.17 0.2 0.14; 0.16 0.17 0.02 0.18; 0 0 0 0.12]
D1=[0 0.88 0.1 0; 0.18 0.07 0.14 0.1; 0.13 0.15 0.15 0.04; 0.31 0.18 0.12 0.27]

disp('Test:');
disp('-----');

disp('[a,A]=MarginalDistributionFromDMAP(D0,D1):');
[a,A]=MarginalDistributionFromDMAP(D0,D1);
disp(a);
disp(A);

assert(length(a)==size(D0,2) && CheckDPHRepresentation(a,A), 'MarginalDistributionFromDMAP returned a wrong DPH representation!');

disp('----------------------------------------------------------------------------');
help MarginalMomentsFromDMAP

disp('Input:');
disp('------');

D0=[0 0.02 0 0; 0 0.17 0.2 0.14; 0.16 0.17 0.02 0.18; 0 0 0 0.12]
D1=[0 0.88 0.1 0; 0.18 0.07 0.14 0.1; 0.13 0.15 0.15 0.04; 0.31 0.18 0.12 0.27]

disp('Test:');
disp('-----');

disp('moms=MarginalMomentsFromDMAP(D0,D1):');
moms=MarginalMomentsFromDMAP(D0,D1);
disp(moms);

assert(length(moms)==2*size(D0,2)-1 && CheckMoments(moms), 'MarginalMomentsFromDMAP returned wrong moments!');

disp('----------------------------------------------------------------------------');
help MarginalDistributionFromDRAP

disp('Input:');
disp('------');

H0=[0 0 0.13; 0 0.6 0.18; 0.31 0.26 0.02]
H1=[0 1 -0.13; 0 0.18 0.04; 0.03 0.09 0.29]

disp('Test:');
disp('-----');

disp('[a,A]=MarginalDistributionFromDRAP(H0,H1):');
[a,A]=MarginalDistributionFromDRAP(H0,H1);
disp(a);
disp(A);

assert(length(a)==size(H0,2) && CheckMGRepresentation(a,A), 'MarginalDistributionFromDRAP returned a wrong MG representation!');

disp('----------------------------------------------------------------------------');
help MarginalMomentsFromDRAP

disp('Input:');
disp('------');

H0=[0 0 0.13; 0 0.6 0.18; 0.31 0.26 0.02]
H1=[0 1 -0.13; 0 0.18 0.04; 0.03 0.09 0.29]

disp('Test:');
disp('-----');

disp('moms=MarginalMomentsFromDRAP(H0,H1):');
moms=MarginalMomentsFromDRAP(H0,H1);
disp(moms);

assert(length(moms)==2*size(H0,2)-1 && CheckMoments(moms), 'MarginalMomentsFromDRAP returned wrong moments!');

disp('----------------------------------------------------------------------------');
help MarginalDistributionFromDMMAP

D0=[0.34 0 0; 0.06 0.05 0.03; 0.11 0.13 0]
D1=[0.3 0 0; 0.16 0.18 0.05; 0.15 0.04 0.09]
D2=[0 0.01 0; 0.1 0.07 0.08; 0.13 0.12 0.13]
D3=[0.35 0 0; 0 0.18 0.04; 0.06 0.03 0.01]

disp('Test:');
disp('-----');

disp('[a,A]=MarginalDistributionFromDMMAP({D0,D1,D2,D3}):');
[a,A]=MarginalDistributionFromDMMAP({D0,D1,D2,D3});
disp(a);
disp(A);

assert(length(a)==size(D0,2) && CheckDPHRepresentation(a,A), 'MarginalDistributionFromDMMAP returned a wrong DPH representation!');

disp('----------------------------------------------------------------------------');
help MarginalMomentsFromDMMAP

D0=[0.34 0 0; 0.06 0.05 0.03; 0.11 0.13 0]
D1=[0.3 0 0; 0.16 0.18 0.05; 0.15 0.04 0.09]
D2=[0 0.01 0; 0.1 0.07 0.08; 0.13 0.12 0.13]
D3=[0.35 0 0; 0 0.18 0.04; 0.06 0.03 0.01]

disp('Test:');
disp('-----');

disp('moms=MarginalMomentsFromDMMAP({D0,D1,D2,D3}):');
moms=MarginalMomentsFromDMMAP({D0,D1,D2,D3});
disp(moms);

assert(length(moms)==2*size(D0,2)-1 && CheckMoments(moms), 'MarginalMomentsFromDMMAP returned wrong moments!');

disp('----------------------------------------------------------------------------');
help MarginalDistributionFromDMRAP

disp('Input:');
disp('------');

H0=[0.15 0.2 0.18; -0.23 0.17 0.22; 0.19 0.15 0.16]
H1=[0.01 0.08 0.16; 0.02 0.2 0.07; 0.02 0.15 0.17]
H2=[0.14 0.07 0.01; 0.19 0.02 0.34; 0.06 0.1 0]

disp('Test:');
disp('-----');

disp('[a,A]=MarginalDistributionFromDMRAP({H0,H1,H2}):');
[a,A]=MarginalDistributionFromDMRAP({H0,H1,H2});
disp(a);
disp(A);

assert(length(a)==size(H0,2) && CheckMGRepresentation(a,A), 'MarginalDistributionFromDMRAP returned a wrong MG representation!');

disp('----------------------------------------------------------------------------');
help MarginalMomentsFromDMRAP

disp('Input:');
disp('------');

H0=[0.15 0.2 0.18; -0.23 0.17 0.22; 0.19 0.15 0.16]
H1=[0.01 0.08 0.16; 0.02 0.2 0.07; 0.02 0.15 0.17]
H2=[0.14 0.07 0.01; 0.19 0.02 0.34; 0.06 0.1 0]

disp('Test:');
disp('-----');

disp('moms=MarginalMomentsFromDMRAP({H0,H1,H2}):');
moms=MarginalMomentsFromDMRAP({H0,H1,H2});
disp(moms);

assert(length(moms)==2*size(H0,2)-1 && CheckMoments(moms), 'MarginalMomentsFromDMRAP returned wrong moments!');

disp('----------------------------------------------------------------------------');
help LagCorrelationsFromDMAP

disp('Input:');
disp('------');

D0=[0 0.02 0 0; 0 0.17 0.2 0.14; 0.16 0.17 0.02 0.18; 0 0 0 0.12]
D1=[0 0.88 0.1 0; 0.18 0.07 0.14 0.1; 0.13 0.15 0.15 0.04; 0.31 0.18 0.12 0.27]

disp('Test:');
disp('-----');

disp('LagCorrelationsFromDMAP(D0,D1,3):');
corr = LagCorrelationsFromDMAP(D0,D1,3);
disp(corr);

assert(length(corr)==3 && all(corr<1) && all(corr>-1), 'LagCorrelationsFromDMAP returned wrong autocorrelation coefficients!');

disp('----------------------------------------------------------------------------');
help LagCorrelationsFromDRAP

disp('Input:');
disp('------');

H0=[0 0 0.13; 0 0.6 0.18; 0.31 0.26 0.02]
H1=[0 1 -0.13; 0 0.18 0.04; 0.03 0.09 0.29]

disp('Test:');
disp('-----');

disp('LagCorrelationsFromDRAP(H0,H1,3):');
corr = LagCorrelationsFromDRAP(H0,H1,3);
disp(corr);

assert(length(corr)==3 && all(corr<1) && all(corr>-1), 'LagCorrelationsFromDRAP returned wrong autocorrelation coefficients!');

disp('----------------------------------------------------------------------------');
help LagkJointMomentsFromDMAP

disp('Input:');
disp('------');

D0=[0 0.02 0 0; 0 0.17 0.2 0.14; 0.16 0.17 0.02 0.18; 0 0 0 0.12]
D1=[0 0.88 0.1 0; 0.18 0.07 0.14 0.1; 0.13 0.15 0.15 0.04; 0.31 0.18 0.12 0.27]

disp('Test:');
disp('-----');

disp('LagkJointMomentsFromDMAP(D0,D1,2,1):');
Nm=LagkJointMomentsFromDMAP(D0,D1,4,1);
disp(Nm);

moms=MarginalMomentsFromDMAP(D0,D1,4);

disp('Correlation from joint moments:');
cjm=zeros(3,1);
for i=1:3
    Nx=LagkJointMomentsFromDMAP(D0,D1,1,i);
    cjm(i) = (Nx(2,2)-moms(1)^2) / (moms(2)-moms(1)^2);
end
disp(cjm)
corr = LagCorrelationsFromDMAP(D0,D1,3);

assert(all(all(Nm>0)) && norm(moms-Nm(1,2:end))<1e-13 && norm(moms'-Nm(2:end,1))<1e-13 && norm(corr-cjm)<1e-13, 'Joint moment matrix is invalid!');

disp('----------------------------------------------------------------------------');
help LagkJointMomentsFromDRAP

disp('Input:');
disp('------');

H0=[0 0 0.13; 0 0.6 0.18; 0.31 0.26 0.02]
H1=[0 1 -0.13; 0 0.18 0.04; 0.03 0.09 0.29]

disp('Test:');
disp('-----');

disp('LagkJointMomentsFromDRAP(H0,H1,2,1):');
Nm=LagkJointMomentsFromDRAP(H0,H1,4,1);
disp(Nm);

moms=MarginalMomentsFromDRAP(H0,H1,4);

disp('Correlation from joint moments:');
cjm=zeros(3,1);
for i=1:3
    Nx=LagkJointMomentsFromDRAP(H0,H1,1,i);
    cjm(i) = (Nx(2,2)-moms(1)^2) / (moms(2)-moms(1)^2);
end
disp(cjm)
corr = LagCorrelationsFromDRAP(H0,H1,3);

assert(all(all(Nm>0)) && norm(moms-Nm(1,2:end))<1e-13 && norm(moms'-Nm(2:end,1))<1e-13 && norm(corr-cjm)<1e-13, 'Joint moment matrix is invalid!');

disp('----------------------------------------------------------------------------');
help LagkJointMomentsFromDMMAP

disp('Input:');
disp('------');

D0=[0.34 0 0; 0.06 0.05 0.03; 0.11 0.13 0]
D1=[0.3 0 0; 0.16 0.18 0.05; 0.15 0.04 0.09]
D2=[0 0.01 0; 0.1 0.07 0.08; 0.13 0.12 0.13]
D3=[0.35 0 0; 0 0.18 0.04; 0.06 0.03 0.01]

disp('Test:');
disp('-----');

disp('LagkJointMomentsFromMMAP({D0,D1,D2,D3},3,1)');
Nm=LagkJointMomentsFromDMMAP({D0,D1,D2,D3},3,1);

disp('Moments of arrival type 1, Nm{1}:');
disp(Nm{1});
disp('Moments of arrival type 2, Nm{2}:');
disp(Nm{2});
disp('Moments of arrival type 2, Nm{3}:');
disp(Nm{3});

assert(length(Nm)==3 && norm(SumMatrixList(Nm)-LagkJointMomentsFromDMAP(D0,D1+D2+D3,3,1))<1e-13,'Joint moment matrix is invalid!');

disp('----------------------------------------------------------------------------');
help LagkJointMomentsFromDMRAP

disp('Input:');
disp('------');

H0=[0.15 0.2 0.18; -0.23 0.17 0.22; 0.19 0.15 0.16]
H1=[0.01 0.08 0.16; 0.02 0.2 0.07; 0.02 0.15 0.17]
H2=[0.14 0.07 0.01; 0.19 0.02 0.34; 0.06 0.1 0]

disp('Test:');
disp('-----');

disp('LagkJointMomentsFromMRAP({H0,H1,H2},3,2)');
Nm=LagkJointMomentsFromDMRAP({H0,H1,H2},3,2);

disp('Moments of arrival type 1, Nm{1}:');
disp(Nm{1});
disp('Moments of arrival type 2, Nm{2}:');
disp(Nm{2});

assert(length(Nm)==2 && norm(SumMatrixList(Nm)-LagkJointMomentsFromDRAP(H0,H1+H2,3,2))<1e-13,'Joint moment matrix is invalid!');

disp('----------------------------------------------------------------------------');
help RandomDMAP

disp('Test:');
disp('-----');

disp('[D0,D1]=RandomDMAP(4,5.62,10):');
[D0,D1]=RandomDMAP(4,5.62,10);
disp(D0);
disp(D1);

disp('Check the mean of the obtained DMAP:');
m = MarginalMomentsFromDMAP(D0,D1,1);
disp(m);

assert(CheckDMAPRepresentation(D0,D1), 'RandomDMAP failed to return a valid DMAP representation!');
assert(max(abs(m-5.62))<1e-14, 'RandomDMAP failed to match the given mean value!');

disp('----------------------------------------------------------------------------');
help RandomDMMAP

disp('Test:');
disp('-----');

disp('D=RandomDMMAP(4,3,5.62,10):');
D=RandomDMMAP(4,3,5.62,10);
disp(D{1});
disp(D{2});
disp(D{3});

disp('Check the mean of the obtained DMMAP:');
m = MarginalMomentsFromDMMAP(D,1);
disp(m);

assert(CheckDMMAPRepresentation(D), 'RandomDMMAP failed to return a valid DMMAP representation!');
assert(max(abs(m-5.62))<1e-14, 'RandomDMMAP failed to match the given mean value!');


disp('----------------------------------------------------------------------------');
help CheckDMAPRepresentation

disp('Input:');
disp('------');

D0=[0 0.02 0; 0 0.17 0.2; 0.16 0.17 0.02]
D1=[0 0.88 0.1 0; 0.18 0.07 0.14 0.1; 0.13 0.15 0.15 0.04; 0.31 0.18 0.12 0.27]

disp('Test:');
disp('-----');

disp('CheckDMAPRepresentation(D0,D1):');
flag=CheckDMAPRepresentation(D0,D1);
disp(flag);

assert(flag==0, 'CheckDMAPRepresentation failed to detect the incompatible shapes of D0 and D1!');

disp('Input:');
disp('------');

D0=[0 0.02 0; 0 0.17 0.2; 0.16 0.17 0.02]
D1=[0 0.88 0.1; 0.18 0.07 0.14; 0.13 0.15 0.15]

disp('Test:');
disp('-----');

disp('CheckDMAPRepresentation(D0,D1):');
flag=CheckDMAPRepresentation(D0,D1);
disp(flag);

assert(flag==0, 'CheckDMAPRepresentation failed to detect invalid rowsums!');

disp('Input:');
disp('------');

D0=[0 0.02 0 0; 0 0.17 0.2 0.14; 0.16 0.17 0.02 0.18; 0 0 0 0.12];
D1=[0 0.88 0.1 0; 0.18 0.07 0.14 0.1; 0.13 0.15 0.15 0.04; 0.31 0.18 0.12 0.27];

disp('Test:');
disp('-----');

disp('CheckDMAPRepresentation(D0,D1):');
flag=CheckDMAPRepresentation(D0,D1);
disp(flag);

assert(flag==1, 'CheckDMAPRepresentation failed to recognize a valid DMAP representation!');

disp('----------------------------------------------------------------------------');
help CheckDRAPRepresentation

disp('Input:');
disp('------');

H0=[0 0 0.13; 0 0.6 0.18; 0.31 0.26 0.02; 0.2 0 0]
H1=[0 1 -0.13; 0 0.18 0.04; 0.03 0.09 0.29; 0 0.8 0]


disp('Test:');
disp('-----');

disp('CheckDRAPRepresentation(H0,H1):')
flag=CheckDRAPRepresentation(H0,H1);
disp(flag);

assert(flag==0, 'CheckDRAPRepresentation failed to detect the incompatible shapes of H0 and H1!');

disp('Input:');
disp('------');

H0=[0.2 0 0.13; 0 0.6 0.18; 0.31 0.26 0.02]
H1=[0 1 -0.13; 0 0.18 0.04; 0.03 0.09 0.29]

disp('Test:');
disp('-----');

disp('CheckDRAPRepresentation(H0,H1):')
flag=CheckDRAPRepresentation(H0,H1);
disp(flag);

assert(flag==0, 'CheckDRAPRepresentation failed to detect invalid rowsums!');

disp('Input:');
disp('------');

x = 15;
H0=[0 0 x; 0 0.6 0.18; 0.31 0.26 0.02]
H1=[0 1 -x; 0 0.18 0.04; 0.03 0.09 0.29]

disp('Test:');
disp('-----');

disp('CheckDRAPRepresentation(H0,H1):')
flag=CheckDRAPRepresentation(H0,H1);
disp(flag);

assert(flag==0, 'CheckDRAPRepresentation failed to detect invalid eigenvalues!');

disp('Input:');
disp('------');

H0 = [0 0.5 0.1; 0 -1.4 3.1; 0.67 0.0 0.4]
H1 = [0 0.4 0; 0 -0.2 -0.5; 0.3 -0.7 0.33]

disp('Test:');
disp('-----');

disp('CheckDRAPRepresentation(H0,H1):')
flag=CheckDRAPRepresentation(H0,H1);
disp(flag);

assert(flag==0, 'CheckDRAPRepresentation failed to detect invalid eigenvalues!');

disp('Input:');
disp('------');

H0=[0 0 0.13; 0 0.6 0.18; 0.31 0.26 0.02]
H1=[0 1 -0.13; 0 0.18 0.04; 0.03 0.09 0.29]

disp('Test:');
disp('-----');

disp('CheckDRAPRepresentation(H0,H1):')
flag=CheckDRAPRepresentation(H0,H1);
disp(flag);

assert(flag==1, 'CheckDRAPRepresentation failed to recognize a valid DRAP representation!');

disp('----------------------------------------------------------------------------');
help CheckDMMAPRepresentation

disp('Input:');
disp('------');

D0=[0.34 0 0; 0.06 0.05 0.03; 0.11 0.13 0]
D1=[0.3 0 0; 0.16 0.18 0.05; 0.15 0.04 0.09]
D2=[0 0.01 0; 0.1 0.07 0.08; 0.13 0.12 0.13]
D3=[0.35 0 0; 0 0.18 0.04; 0.06 0.03 0.01]

disp('Test:');
disp('-----');

disp('CheckDMMAPRepresentation({D0,D1,D2,D3}):');
flag=CheckDMMAPRepresentation({D0,D1,D2,D3});
disp(flag);

assert(flag==1, 'CheckDMMAPRepresentation failed to recognize a valid DMMAP representation!');

disp('----------------------------------------------------------------------------');
help CheckDMRAPRepresentation

disp('Input:');
disp('------');

H0=[0.15 0.2 0.18; -0.23 0.17 0.22; 0.19 0.15 0.16]
H1=[0.01 0.08 0.16; 0.02 0.2 0.07; 0.02 0.15 0.17]
H2=[0.14 0.07 0.01; 0.19 0.02 0.34; 0.06 0.1 0]

disp('Test:');
disp('-----');

disp('CheckDMRAPRepresentation({H0,H1,H2}):');
flag=CheckDMRAPRepresentation({H0,H1,H2});
disp(flag);

assert(flag==1, 'CheckDMRAPRepresentation failed to recognize a valid DMRAP representation!');

disp('----------------------------------------------------------------------------');
help DRAPFromMoments

disp('Input:');
disp('------');

G0=[0 0.02 0; 0 0.17 0.2; 0.16 0.17 0.24];
G1=[0 0.88 0.1; 0.42 0.07 0.14; 0.13 0.15 0.15];

moms=MarginalMomentsFromDRAP(G0,G1,5)
Nm=LagkJointMomentsFromDRAP(G0,G1,2,1)

disp('Test:');
disp('-----');

disp('[H0,H1]=DRAPFromMoments(moms,Nm):');
[H0,H1]=DRAPFromMoments(moms,Nm);
disp(H0);
disp(H1);

rmoms=MarginalMomentsFromDRAP(H0,H1,5,1e-12)
rNm=LagkJointMomentsFromDRAP(H0,H1,2,1,1e-12)

assert(norm(moms-rmoms)<1e-11 && norm(Nm-rNm)<1e-12, 'The moments and joint moments returned by DRAPFromMoments are not the same as given!');

disp('----------------------------------------------------------------------------');
help DMRAPFromMoments

disp('Input:');
disp('------');

G0=[0.34 0 0; 0.06 0.05 0.03; 0.11 0.13 0]
G1=[0.3 0 0; 0.16 0.18 0.05; 0.15 0.04 0.09]
G2=[0 0.01 0; 0.1 0.07 0.08; 0.13 0.12 0.13]
G3=[0.35 0 0; 0 0.18 0.04; 0.06 0.03 0.01]
G={G0,G1,G2,G3};

moms=MarginalMomentsFromDMRAP(G,5)
Nm=LagkJointMomentsFromDMRAP(G,2,1);
Nm1=Nm{1}
Nm2=Nm{2}
Nm3=Nm{3}

disp('Test:');
disp('-----');

disp('H=DMRAPFromMoments(moms,Nm):');
H=DMRAPFromMoments(moms,Nm);
disp(H{1});
disp(H{2});
disp(H{3});
disp(H{4});

rmoms=MarginalMomentsFromDMRAP(H,5,1e-11)
rNm=LagkJointMomentsFromDMRAP(H,2,1,1e-11);
rNm1=rNm{1}
rNm2=rNm{2}
rNm3=rNm{3}

assert(norm(moms-rmoms)<1e-10 && norm(Nm1-rNm1)<1e-10 && norm(Nm2-rNm2)<1e-10 && norm(Nm3-rNm3)<1e-10, 'The moments and joint moments returned by DMRAPFromMoments are not the same as given!');

disp('----------------------------------------------------------------------------');
help DMAPFromDRAP

disp('Input:');
disp('------');

x=0.13;
H0=[0 0 x; 0 0.6 0.18; 0.31 0.26 0.02]
H1=[0 1 -x; 0 0.18 0.04; 0.03 0.09 0.29]

disp('Test:');
disp('-----');

disp('[D0,D1]=DMAPFromDRAP(H0,H1):');
[D0,D1]=DMAPFromDRAP(H0,H1);
disp(D0);
disp(D1);

err = norm(LagkJointMomentsFromDRAP(H0,H1,3,1,1e-12)-LagkJointMomentsFromDRAP(D0,D1,3,1,1e-12));
assert(err<1e-10, 'The DRAP returned by DMAPFromDRAP is not similar to the input!');
assert(CheckDMAPRepresentation(D0,D1), 'The result of DMAPFromDRAP is not a DMAP, as it should be!');

disp('----------------------------------------------------------------------------');
help DMMAPFromDMRAP

disp('Input:');
disp('------');

H0=[0.15 0.2 0.18; -0.20 0.17 0.22; 0.19 0.15 0.16]
H1=[0.01 0.08 0.16; 0.02 0.2 0.07; 0.02 0.15 0.17]
H2=[0.14 0.07 0.01; 0.19 0.02 0.31; 0.06 0.1 0]
H = {H0,H1,H2};

moms=MarginalMomentsFromDMRAP(H);
jmom=LagkJointMomentsFromDMRAP(H,3,1);

disp('Test:');
disp('-----');

disp('D=DMMAPFromDMRAP({H0,H1,H2}):');
D=DMMAPFromDMRAP(H);
disp(D{1})
disp(D{2})
disp(D{3})

rmoms=MarginalMomentsFromDMMAP(D);
rjmom=LagkJointMomentsFromDMMAP(D,3,1);

err = norm(rjmom{1}-jmom{1}) + norm(rjmom{2}-jmom{2});
assert(err<1e-12, 'The DMMAP returned by DMMAPFromDMRAP is not similar to the input!');
assert(CheckDMMAPRepresentation(D), 'The result of DMMAPFromDMRAP is not a DMMAP, as it should be!');

disp('----------------------------------------------------------------------------');
help CanonicalFromDMAP2

disp('Input:');
disp('------');

D0=[0.46 0.28; 0.35 0.23]
D1=[0.08 0.18; 0.14 0.28]

disp('Test:');
disp('-----');

disp('[H0,H1]=CanonicalFromDMAP2(D0,D1):');
[H0,H1]=CanonicalFromDMAP2(D0,D1);
disp(H0);
disp(H1);

C=SimilarityMatrix(H0,D0);
err = norm(H0*C-C*D0) + norm(H1*C-C*D1);

assert(CheckDMAPRepresentation(H0,H1), 'The result of CanonicalFromDMAP2 is not a valid DMAP representation!');
assert(err<1e-12, 'The DMAP returned by CanonicalFromDMAP2 is not similar to the input!');

disp('Input:');
disp('------');

D0=[0.26 0.28; 0.35 0.23]
D1=[0.28 0.18; 0.14 0.28]

disp('Test:');
disp('-----');

disp('[H0,H1]=CanonicalFromDMAP2(D0,D1):');
[H0,H1]=CanonicalFromDMAP2(D0,D1);
disp(H0);
disp(H1);

C=SimilarityMatrix(H0,D0);
err = norm(H0*C-C*D0) + norm(H1*C-C*D1);

assert(CheckDMAPRepresentation(H0,H1), 'The result of CanonicalFromDMAP2 is not a valid DMAP representation!');
assert(err<1e-12, 'The DMAP returned by CanonicalFromDMAP2 is not similar to the input!');

disp('Input:');
disp('------');

D0=[0.14 0.34; 0.35 0.23]
D1=[0.22 0.3; 0.28 0.14]

disp('Test:');
disp('-----');

disp('[H0,H1]=CanonicalFromDMAP2(D0,D1):');
[H0,H1]=CanonicalFromDMAP2(D0,D1);
disp(H0);
disp(H1);

C=SimilarityMatrix(H0,D0);
err = norm(H0*C-C*D0) + norm(H1*C-C*D1);

assert(CheckDMAPRepresentation(H0,H1), 'The result of CanonicalFromDMAP2 is not a valid DMAP representation!');
assert(err<1e-12, 'The DMAP returned by CanonicalFromDMAP2 is not similar to the input!');

disp('----------------------------------------------------------------------------');
help DMAP2FromMoments

disp('Input:');
disp('------');

D0=[0.2 0.7; 0.6 0.1];
D1=[0.09 0.01; 0.2 0.1];
moms=MarginalMomentsFromDMAP(D0, D1, 3)
corr=LagCorrelationsFromDMAP(D0, D1, 1)

disp('Test:');
disp('-----');

disp('[D0,D1]=DMAP2FromMoments(moms,corr):');
[D0,D1]=DMAP2FromMoments(moms,corr);
disp(D0);
disp(D1);

rmoms=MarginalMomentsFromDMAP(D0, D1, 3)
rcorr=LagCorrelationsFromDMAP(D0, D1, 1)

assert(CheckDMAPRepresentation(D0,D1), 'DMAP2FromMoments returned an invalid DMAP representation!');
assert(norm(moms-rmoms)<1e-11 && norm(corr-rcorr)<1e-11, 'The moments and the correlation returned by DMAP2FromMoments are not the same as given!');

disp('----------------------------------------------------------------------------');
help SamplesFromDMAP

disp('Input:');
disp('------');

D0=[0 0.02 0 0; 0 0.17 0.2 0.14; 0.16 0.17 0.02 0.18; 0 0 0 0.12]
D1=[0 0.88 0.1 0; 0.18 0.07 0.14 0.1; 0.13 0.15 0.15 0.04; 0.31 0.18 0.12 0.27]

disp('Test:');
disp('-----');

disp('x=SamplesFromDMAP(D0,D1,10000)');
x=SamplesFromDMAP(D0,D1,10000);

disp('Moments from the samples:');
mt = MarginalMomentsFromTrace(x,3);
disp(mt);

disp('Moments from the DMAP:');
mm = MarginalMomentsFromDMAP(D0,D1,3);
disp(mm);

disp('----------------------------------------------------------------------------');
help SamplesFromDMMAP

disp('Input:');
disp('------');

D0=[0.34 0 0; 0.06 0.05 0.03; 0.11 0.13 0]
D1=[0.3 0 0; 0.16 0.18 0.05; 0.15 0.04 0.09]
D2=[0 0.01 0; 0.1 0.07 0.08; 0.13 0.12 0.13]
D3=[0.35 0 0; 0 0.18 0.04; 0.06 0.03 0.01]
D = {D0,D1,D2,D3};

disp('Test:');
disp('-----');

disp('x=SamplesFromDMMAP(D,10000)');
x=SamplesFromDMMAP(D,10000);

disp('Moments from the samples:');
mt = MarginalMomentsFromTrace(x(:,1),3);
disp(mt);

disp('Moments from the DMMAP:');
mm = MarginalMomentsFromDMMAP(D,3);
disp(mm);

