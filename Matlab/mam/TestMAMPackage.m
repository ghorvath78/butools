format short g;

disp('---BuTools: test file for Matrix analytic methods ---');

disp('Enable the verbose messages with the BuToolsVerbose flag')
global BuToolsVerbose;
BuToolsVerbose = true

disp('Enable input parameter checking with the BuToolsCheckInput flag')
global BuToolsCheckInput;
BuToolsCheckInput = true

B = [0,0;3,4];
L = [-6,5;3,-12];
F = [1,0;2,0];
L0 = [-6,5;6,-8];

disp('----------------------------------------------------------------------------');
help QBDFundamentalMatrices

disp('Input:');
disp('------');
B
L
F
L0

disp('Test:');
disp('-----');

disp('M = QBDFundamentalMatrices (B,L,F,''RGU''):');
M = QBDFundamentalMatrices (B,L,F,'RGU');
R=M{1}
G=M{2}
U=M{3}

assert(CheckGenerator(U,1)==1, 'QBDFundamentalMatrices: matrix U is not a transient generator!');
assert(all(abs(eig(R))<1), 'QBDFundamentalMatrices: the eigenvalues of matrix R are not inside the unit disc!');
assert(CheckProbMatrix(G)==1, 'QBDFundamentalMatrices: matrix G is not a transition prob. matrix!');

disp('----------------------------------------------------------------------------');
help QBDSolve

disp('Test:');
disp('-----');

disp('[pi0, R] = QBDSolve (B, L, F, L0):');
[pi0, R] = QBDSolve (B, L, F, L0)

assert(sum(pi0)>0 && sum(pi0)<=1 && all(pi0>=0), 'QBDSolve: wrong pi0 vector!');
assert(all(abs(eig(R))<1), 'QBDSolve: the eigenvalues of matrix R are not inside the unit disc!');


disp('----------------------------------------------------------------------------');
help QBDStationaryDistr

disp('Test:');
disp('-----');

disp('pi = QBDStationaryDistr (pi0, R, 5):');
pi = QBDStationaryDistr (pi0, R, 5)

assert(sum(pi)>0 && sum(pi)<=1 && all(pi>=0), 'QBDStationaryDistr: wrong pi vector!');

%==================== M/G/1 type queues ==============================
B0 = [0.1, 0.5; 0.3, 0.4];
B1 = [0, 0.1; 0, 0];
B2 = [0.2, 0; 0, 0.2];
B3 = [0, 0.1; 0.1, 0];
A0 = [0.4, 0.2; 0.3, 0.4];
A1 = [0, 0.1; 0, 0];
A2 = [0, 0.2; 0, 0.2];
A3 = [0.1, 0; 0.1, 0];

B = [B0,B1,B2,B3];
A = [A0,A1,A2,A3];

G = MG1FundamentalMatrix (A)
pi = MG1StationaryDistr (A,B,G,300);

assert(CheckProbMatrix(G)==1, 'MG1FundamentalMatrix: matrix G is not a transition prob. matrix!');
assert(sum(pi)>0 && sum(pi)<=1 && all(pi>=0), 'MG1StationaryDistr: wrong pi vector!');

%==================== G/M/1 type queues ==============================
B0 = [0.7, 0.2; 0.3, 0.6];
B1 = [0.3, 0.4; 0.5, 0.2];
B2 = [0.2, 0.4; 0.1, 0.6];
B3 = [0, 0.1; 0.2, 0];
A0 = [0.1, 0; 0, 0.1];
A1 = [0, 0.2; 0, 0.2];
A2 = [0, 0.1; 0, 0];
A3 = [0.3, 0.2; 0.3, 0.2];
A4 = B3;

B = [B0;B1;B2;B3];
A = [A0,A1,A2,A3,A4];

R = GM1FundamentalMatrix (A)
pi = GM1StationaryDistr (B,R,300);

assert(all(abs(eig(R))<1), 'GM1FundamentalMatrix: the eigenvalues of matrix R are not inside the unit disc!');
assert(sum(pi)>0 && sum(pi)<=1 && all(pi>=0), 'GM1StationaryDistr: wrong pi vector!');

%==================== fluid tests ==============================

Fpp=[-5 1; 2 -3];
Fpm=[2 1 1; 1 0 0];
Fmm=[-8 4 1; 2 -12 3; 2 0 -2];
Fmp=[3 0; 2 5; 0 0];

QA = [-2 2; 5 -5];
RA = diag([3 7]);
NA = size(QA,1);
IA = eye(NA);
QS = [-4 1 3; 6 -8 2; 3 7 -10];
RS = diag([1 7 15]);
NS = size(QS,1);
IS = eye(NS);

Q = kron(QA, IS) + kron(IA, QS);
R = kron(RA, IS) - kron(IA, RS);

disp('----------------------------------------------------------------------------');
help FluidFundamentalMatrices

disp('Input:');
disp('------');
Fpp
Fpm
Fmp
Fmm

disp('Test:');
disp('-----');

disp('M = FluidFundamentalMatrices (Fpp, Fpm, Fmp, Fmm, ''PKU''):');
M = FluidFundamentalMatrices (Fpp, Fpm, Fmp, Fmm, 'PKU');
Psi=M{1}
K=M{2}
U=M{3}
assert(CheckGenerator(U)==1, 'FluidFundamentalMatrices: matrix U is not a generator!');
assert(all(eig(K)<0), 'FluidFundamentalMatrices: the eigenvalues of matrix K are not negative!');
assert(all(all(Psi>=0)) && all(all(Psi<=1)) && norm(sum(Psi,2)-1)<1e-14, 'FluidFundamentalMatrices: matrix Psi is not a transition prob. matrix!');

disp('----------------------------------------------------------------------------');
help FluidSolve

disp('Test:');
disp('-----');

disp('[mass0, ini, K, clo] = FluidSolve (Fpp, Fpm, Fmp, Fmm):');
[mass0, ini, K, clo] = FluidSolve (Fpp, Fpm, Fmp, Fmm)

one = sum(mass0) + sum(ini*inv(-K)*clo);
assert(abs(one-1)<1e-14, 'FluidSolve: the integral of the fluid distribution in not one!');

disp('----------------------------------------------------------------------------');
help GeneralFluidSolve

disp('Test:');
disp('-----');

disp('[mass0, ini, K, clo] = GeneralFluidSolve (Q, R):');
[mass0, ini, K, clo] = GeneralFluidSolve (Q, R)

one = sum(mass0) + sum(ini*inv(-K)*clo);
assert(abs(one-1)<1e-14, 'GeneralFluidSolve: the integral of the fluid distribution in not one!');

disp('----------------------------------------------------------------------------');
help FluidStationaryDistr

disp('Test:');
disp('-----');

disp('y = FluidStationaryDistr (mass0, ini, K, clo, (0:0.1:1)''):');
y = FluidStationaryDistr (mass0, ini, K, clo, (0:1:30)')

pi=CTMCSolve(Q)
assert(norm(y(end,:)-pi)<1e-5, 'FluidStationaryDistr: stationary distribution does not converge to the steady state distribution of the phases!');


