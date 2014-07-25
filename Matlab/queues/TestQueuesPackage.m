format short g;

disp('---BuTools: test file for Matrix analytic methods ---');

disp('Enable the verbose messages with the BuToolsVerbose flag')
global BuToolsVerbose;
BuToolsVerbose = true

disp('Enable input parameter checking with the BuToolsCheckInput flag')
global BuToolsCheckInput;
BuToolsCheckInput = true

% ============================= QBD tests ===============================

disp('----------------------------------------------------------------------------');
help QBDQueueQLD
help QBDQueueSTD

disp('Input:');
disp('------');

% QBD example 1
B=[6, 1, 0; 0, 4, 1; 2, 0, 0]
F=[0, 1, 1; 5, 0, 0; 1, 3, 0]
L=[-14, 3, 2; 0, -14, 4; 3, 1, -10]
L0 = L+B

[pi0, R] = QBDSolve (B, L, F, L0);
lambda=sum(pi0*inv(eye(size(R))-R)*F);

disp('Test:');
disp('-----');

disp('[alphap,Ap] = QBDQueueQLD(B, L, F, L0, true):');
[alphap,Ap] = QBDQueueQLD(B, L, F, L0, true)
disp('[alpha,A] = QBDQueueQLD(B, L, F, L0):');
[alpha,A] = QBDQueueQLD(B, L, F, L0)
disp('[betap, Bp] = QBDQueueSTD(B, L, F, L0, true):');
[betap, Bp] = QBDQueueSTD(B, L, F, L0, true)
disp('[beta, B] = QBDQueueSTD(B, L, F, L0):');
[beta, B] = QBDQueueSTD(B, L, F, L0)

assert(CheckMGRepresentation(alpha,A), 'QBDQueueQLD: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'QBDQueueSTD: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'QBDQueueQLD: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'QBDQueueSTD: invalid PH representation of the sojourn time!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
assert(abs(mql-mst*lambda)<1e-12, 'QBDQueue: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'QBDQueueQLD: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'QBDQueueSTD: the ME and PH representations are not equal!');

disp('Input:');
disp('------');

% QBD example 2
B=[6, 1, 0; 0, 4, 1; 2, 0, 0]
F=[0, 0, 0; 5, 0, 0; 1, 3, 0]
L=[-12, 3, 2; 0, -14, 4; 3, 1, -10]
L0 = L+B

[pi0, R] = QBDSolve (B, L, F, L0);
lambda=sum(pi0*inv(eye(size(R))-R)*F);

disp('Test:');
disp('-----');

disp('[alphap,Ap] = QBDQueueQLD(B, L, F, L0, true):');
[alphap,Ap] = QBDQueueQLD(B, L, F, L0, true)
disp('[alpha,A] = QBDQueueQLD(B, L, F, L0):');
[alpha,A] = QBDQueueQLD(B, L, F, L0)
disp('[betap, Bp] = QBDQueueSTD(B, L, F, L0, true):');
[betap, Bp] = QBDQueueSTD(B, L, F, L0, true)
disp('[beta, B] = QBDQueueSTD(B, L, F, L0):');
[beta, B] = QBDQueueSTD(B, L, F, L0)

assert(CheckMGRepresentation(alpha,A), 'QBDQueueQLD: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'QBDQueueSTD: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'QBDQueueQLD: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'QBDQueueSTD: invalid PH representation of the sojourn time!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
assert(abs(mql-mst*lambda)<1e-12, 'QBDQueue: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'QBDQueueQLD: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'QBDQueueSTD: the ME and PH representations are not equal!');

disp('Input:');
disp('------');

% QBD example 3
B=[6, 1, 0; 0, 5, 0; 0, 0, 0]
F=[0, 3, 1; 0, 5, 0; 0, 0, 0]
L=[-16, 3, 2; 0, -14, 4; 3, 1, -4]
L0=[-14, 10, 0; 5, -10, 0; 3, 1, -4]

[pi0, R] = QBDSolve (B, L, F, L0);
lambda=sum(pi0*inv(eye(size(R))-R)*F);

disp('Test:');
disp('-----');

disp('[alphap,Ap] = QBDQueueQLD(B, L, F, L0, true):');
[alphap,Ap] = QBDQueueQLD(B, L, F, L0, true)
disp('[alpha,A] = QBDQueueQLD(B, L, F, L0):');
[alpha,A] = QBDQueueQLD(B, L, F, L0)
disp('[betap, Bp] = QBDQueueSTD(B, L, F, L0, true):');
[betap, Bp] = QBDQueueSTD(B, L, F, L0, true)
disp('[beta, B] = QBDQueueSTD(B, L, F, L0):');
[beta, B] = QBDQueueSTD(B, L, F, L0)

assert(CheckMGRepresentation(alpha,A), 'QBDQueueQLD: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'QBDQueueSTD: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'QBDQueueQLD: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'QBDQueueSTD: invalid PH representation of the sojourn time!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
assert(abs(mql-mst*lambda)<1e-12, 'QBDQueue: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'QBDQueueQLD: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'QBDQueueSTD: the ME and PH representations are not equal!');

disp('Input:');
disp('------');

% QBD example 4
B = [0,0; 3,4]
L = [-6,5; 3,-12]
F = [1,0; 2,0]
L0 = [-6,5; 6,-8]

[pi0, R] = QBDSolve (B, L, F, L0);
lambda=sum(pi0*inv(eye(size(R))-R)*F);

disp('Test:');
disp('-----');

disp('[alphap,Ap] = QBDQueueQLD(B, L, F, L0, true):');
[alphap,Ap] = QBDQueueQLD(B, L, F, L0, true)
disp('[alpha,A] = QBDQueueQLD(B, L, F, L0):');
[alpha,A] = QBDQueueQLD(B, L, F, L0)
disp('[betap, Bp] = QBDQueueSTD(B, L, F, L0, true):');
[betap, Bp] = QBDQueueSTD(B, L, F, L0, true)
disp('[beta, B] = QBDQueueSTD(B, L, F, L0):');
[beta, B] = QBDQueueSTD(B, L, F, L0)

assert(CheckMGRepresentation(alpha,A), 'QBDQueueQLD: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'QBDQueueSTD: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'QBDQueueQLD: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'QBDQueueSTD: invalid PH representation of the sojourn time!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
assert(abs(mql-mst*lambda)<1e-12, 'QBDQueue: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'QBDQueueQLD: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'QBDQueueSTD: the ME and PH representations are not equal!');

% ============================= MAP/MAP/1 tests ===============================

disp('----------------------------------------------------------------------------');
help MAPMAP1QLD
help MAPMAP1STD

disp('Input:');
disp('------');

D0=[-8, 2; 1, -3]
D1=[1, 5; 0, 2]

S0=[-10, 4; 0, -7]
S1=[5, 1; 4, 3]

disp('Test:');
disp('-----');

disp('[alphap,Ap] = MAPMAP1QLD(D0,D1,S0,S1,true):');
[alphap,Ap] = MAPMAP1QLD(D0,D1,S0,S1,true)
disp('[alpha,A] = MAPMAP1QLD(D0,D1,S0,S1):');
[alpha,A] = MAPMAP1QLD(D0,D1,S0,S1)
disp('[betap,Bp] = MAPMAP1STD(D0,D1,S0,S1,true):');
[betap,Bp] = MAPMAP1STD(D0,D1,S0,S1,true)
disp('[beta,B] = MAPMAP1STD(D0,D1,S0,S1):');
[beta,B] = MAPMAP1STD(D0,D1,S0,S1)

assert(CheckMGRepresentation(alpha,A), 'MAPMAP1QLD: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'MAPMAP1STD: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'MAPMAP1QLD: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'MAPMAP1STD: invalid PH representation of the sojourn time!');

% cross-check
IA = eye(size(D0,1));
IS = eye(size(S0,1));
[gamma, G] = QBDQueueSTD (kron(IA,S1), kron(D0,IS)+kron(IA,S0), kron(D1,IS), kron(D0,IS));
msmall = MomentsFromME(beta,B,5);
mlarge = MomentsFromME(gamma,G,5);
assert(norm((msmall-mlarge)./msmall)<1e-12, 'MAPMAP1STD: Large and small model does not give the same results!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
lambda = 1/MarginalMomentsFromMAP(D0,D1,1);
assert(abs(mql-mst*lambda)<1e-12, 'MAPMAP1: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'MAPMAP1QLD: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'MAPMAP1STD: the ME and PH representations are not equal!');

disp('Input:');
disp('------');

D0 = [-8, 1, 2; 0, -6, 4; 3 0 -3]
D1 = [4, 1, 0; 0, 2, 0; 0, 0, 0]

disp('Test:');
disp('-----');

disp('[alphap,Ap] = MAPMAP1QLD(D0,D1,S0,S1,true):');
[alphap,Ap] = MAPMAP1QLD(D0,D1,S0,S1,true)
disp('[alpha,A] = MAPMAP1QLD(D0,D1,S0,S1):');
[alpha,A] = MAPMAP1QLD(D0,D1,S0,S1)
disp('[betap,Bp] = MAPMAP1STD(D0,D1,S0,S1,true):');
[betap,Bp] = MAPMAP1STD(D0,D1,S0,S1,true)
disp('[beta,B] = MAPMAP1STD(D0,D1,S0,S1):');
[beta,B] = MAPMAP1STD(D0,D1,S0,S1)

assert(CheckMGRepresentation(alpha,A), 'MAPMAP1QLD: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'MAPMAP1STD: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'MAPMAP1QLD: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'MAPMAP1STD: invalid PH representation of the sojourn time!');

% cross-check
IA = eye(size(D0,1));
IS = eye(size(S0,1));
[gamma, G] = QBDQueueSTD (kron(IA,S1), kron(D0,IS)+kron(IA,S0), kron(D1,IS), kron(D0,IS));
msmall = MomentsFromME(beta,B,5);
mlarge = MomentsFromME(gamma,G,5);
assert(norm((msmall-mlarge)./msmall)<1e-12, 'MAPMAP1STD: Large and small model does not give the same results!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
lambda = 1/MarginalMomentsFromMAP(D0,D1,1);
assert(abs(mql-mst*lambda)<1e-12, 'MAPMAP1: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'MAPMAP1QLD: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'MAPMAP1STD: the ME and PH representations are not equal!');

disp('Input:');
disp('------');

S0 = [-10, 4, 0; 5, -7, 2; 1, 2, -8]
S1 = [0, 0, 6; 0, 0, 0; 0, 3, 2]

D0 = [-8, 1, 2; 0, -6, 4; 3 0 -3]
D1 = [4, 1, 0; 0, 0, 2; 0, 0, 0]

disp('Test:');
disp('-----');

disp('[alphap,Ap] = MAPMAP1QLD(D0,D1,S0,S1,true):');
[alphap,Ap] = MAPMAP1QLD(D0,D1,S0,S1,true)
disp('[alpha,A] = MAPMAP1QLD(D0,D1,S0,S1):');
[alpha,A] = MAPMAP1QLD(D0,D1,S0,S1)
disp('[betap,Bp] = MAPMAP1STD(D0,D1,S0,S1,true):');
[betap,Bp] = MAPMAP1STD(D0,D1,S0,S1,true)
disp('[beta,B] = MAPMAP1STD(D0,D1,S0,S1):');
[beta,B] = MAPMAP1STD(D0,D1,S0,S1)

assert(CheckMGRepresentation(alpha,A), 'MAPMAP1QLD: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'MAPMAP1STD: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'MAPMAP1QLD: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'MAPMAP1STD: invalid PH representation of the sojourn time!');

% cross-check
IA = eye(size(D0,1));
IS = eye(size(S0,1));
[gamma, G] = QBDQueueSTD (kron(IA,S1), kron(D0,IS)+kron(IA,S0), kron(D1,IS), kron(D0,IS));
msmall = MomentsFromME(beta,B,5);
mlarge = MomentsFromME(gamma,G,5);
assert(norm((msmall-mlarge)./msmall)<1e-12, 'MAPMAP1STD: Large and small model does not give the same results!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
lambda = 1/MarginalMomentsFromMAP(D0,D1,1);
assert(abs(mql-mst*lambda)<1e-12, 'MAPMAP1: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'MAPMAP1QLD: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'MAPMAP1STD: the ME and PH representations are not equal!');

% ============================= MAP/PH/1 tests ===============================

disp('----------------------------------------------------------------------------');
help MAPPH1QLD
help MAPPH1STD

disp('Input:');
disp('------');

D0 = [-8, 1, 2; 0, -6, 4; 3 0 -3]
D1 = [4, 1, 0; 0, 0, 2; 0, 0, 0]
sigma = [0.2, 0.7, 0.1]
S = [-10, 4, 0; 5, -7, 2; 1, 2, -8]

disp('Test:');
disp('-----');

disp('[alphap,Ap] = MAPPH1QLD(D0,D1,sigma,S,true):');
[alphap,Ap] = MAPPH1QLD(D0,D1,sigma,S,true)
disp('[alpha,A] = MAPPH1QLD(D0,D1,sigma,S):');
[alpha,A] = MAPPH1QLD(D0,D1,sigma,S)
disp('[betap,Bp] = MAPPH1STD(D0,D1,sigma,S,true):');
[betap,Bp] = MAPPH1STD(D0,D1,sigma,S,true)
disp('[beta,B] = MAPPH1STD(D0,D1,sigma,S):');
[beta,B] = MAPPH1STD(D0,D1,sigma,S)

assert(CheckMGRepresentation(alpha,A), 'MAPPH1QLD: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'MAPPH1STD: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'MAPPH1QLD: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'MAPPH1STD: invalid PH representation of the sojourn time!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
lambda = 1/MarginalMomentsFromMAP(D0,D1,1);
assert(abs(mql-mst*lambda)<1e-12, 'MAPPH1: Little formula does not hold!');

% ============================= PH/PH/1 tests ===============================

disp('----------------------------------------------------------------------------');
help PHPH1QLD
help PHPH1STD

disp('Input:');
disp('------');

delta = [0.5, 0.1, 0.4]
D = [-8, 1, 2; 0, -6, 4; 3 0 -3]
sigma = [0.2, 0.7, 0.1]
S = [-10, 4, 0; 5, -7, 2; 1, 2, -8]

disp('Test:');
disp('-----');

disp('[alphap,Ap] = PHPH1QLD(delta,D,sigma,S,true):');
[alphap,Ap] = PHPH1QLD(delta,D,sigma,S,true)
disp('[alpha,A] = PHPH1QLD(delta,D,sigma,S):');
[alpha,A] = PHPH1QLD(delta,D,sigma,S)
disp('[betap,Bp] = PHPH1STD(delta,D,sigma,S,true):');
[betap,Bp] = PHPH1STD(delta,D,sigma,S,true)
disp('[beta,B] = PHPH1STD(delta,D,sigma,S):');
[beta,B] = PHPH1STD(delta,D,sigma,S)

assert(CheckMGRepresentation(alpha,A), 'PHPH1QLD: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'PHPH1STD: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'PHPH1QLD: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'PHPH1STD: invalid PH representation of the sojourn time!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
lambda = 1/MomentsFromPH(delta,D,1);
assert(abs(mql-mst*lambda)<1e-12, 'PHPH1: Little formula does not hold!');

% ============================= Fluid tests ===============================

disp('----------------------------------------------------------------------------');
help FluidQueueQLD
help FluidQueueSTD

disp('Input:');
disp('------');
Q = [-9 2 4 0 1 2; 6 -25 5 3 7 4; 1 3 -4 0 0 0; 0 0 0 -8 3 5; 7 3 0 2 -13 1; 7 8 0 3 8 -26]
Rin = diag([4 2 1 0 0 3])
Rout = diag([6 2 0 0 3 2])
lambda = sum(CTMCSolve(Q)*Rin);

disp('Test:');
disp('-----');

disp('[alpha, A] = FluidQueueQLD(Q, Rin, Rout):');
[alpha, A] = FluidQueueQLD(Q, Rin, Rout)
disp('[alphap, Ap] = FluidQueueQLD(Q, Rin, Rout, true):');
[alphap, Ap] = FluidQueueQLD(Q, Rin, Rout, [], true)
disp('[beta, B] = FluidQueueSTD(Q, Rin, Rout):');
[beta, B] = FluidQueueSTD(Q, Rin, Rout);
disp('[betap, Bp] = FluidQueueSTD(Q, Rin, Rout, true):');
[betap, Bp] = FluidQueueSTD(Q, Rin, Rout, [], true);

assert(CheckMERepresentation(alpha,A), 'FluidQueueQLD: invalid ME representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'FluidQueueSTD: invalid ME representation of the sojourn time!');
assert(CheckPHRepresentation(alphap,Ap), 'FluidQueueQLD: invalid PH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'FluidQueueSTD: invalid PH representation of the sojourn time!');

% check Little formula
mql = MomentsFromME(alpha,A,1);
mst = MomentsFromME(beta,B,1);
assert(abs(mql-mst*lambda)<1e-12, 'FluidQueueSTD: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm(MomentsFromPH(alphap,Ap,5)-MomentsFromME(alpha,A,5))<1e-12, 'FluidQueueQLD: the ME and PH representations are not equal!');
assert(norm(MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))<1e-12, 'FluidQueueSTD: the ME and PH representations are not equal!');

% ============================= FluFlu tests ===============================

disp('----------------------------------------------------------------------------');
help FluFluQLD
help FluFluSTD

disp('Input:');
disp('------');

Qin = [-2 1 1; 2 -5 3; 4 0 -4]
Rin = diag([3 7 0])
lambda = sum(CTMCSolve(Qin)*Rin);

Qout = [-4 1 3; 6 -8 2; 3 7 -10]
Rout = diag([1 7 15])

disp('Test:');
disp('-----');

disp('[alphap, Ap] = FluFluQLD(Qin, Rin, Qout, Rout, false, true):');
[alphap, Ap] = FluFluQLD(Qin, Rin, Qout, Rout, false, true)
disp('[alpha, A] = FluFluQLD(Qin, Rin, Qout, Rout, false):');
[alpha, A] = FluFluQLD(Qin, Rin, Qout, Rout, false)
disp('[betap, Bp] = FluFluSTD(Qin, Rin, Qout, Rout, false, true):');
[betap, Bp] = FluFluSTD(Qin, Rin, Qout, Rout, false, true)
disp('[beta, B] = FluFluSTD(Qin, Rin, Qout, Rout, false):');
[beta, B] = FluFluSTD(Qin, Rin, Qout, Rout, false)

assert(CheckMERepresentation(alpha,A), 'FluFluQLD: invalid ME representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'FluFluSTD: invalid ME representation of the sojourn time!');
assert(CheckPHRepresentation(alphap,Ap), 'FluFluQLD: invalid PH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'FluFluSTD: invalid PH representation of the sojourn time!');

% cross-check
Iin = eye(size(Qin));
Iout = eye(size(Qout));
[gamma, G] = FluidQueueSTD (kron(Qin,Iout)+kron(Iin,Qout), kron(Rin,Iout), kron(Iin,Rout));
msmall = MomentsFromME(beta,B,5);
mlarge = MomentsFromME(gamma,G,5);
assert(norm((msmall-mlarge)./msmall)<1e-12, 'FluFluSTD: Large and small model does not give the same results!');

% Little formula
mql = MomentsFromME(alpha,A,1);
mst = MomentsFromME(beta,B,1);
assert(abs(mql-mst*lambda)<1e-12, 'FluFluSTD: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm(MomentsFromPH(alphap,Ap,5)-MomentsFromME(alpha,A,5))<1e-12, 'FluFluQLD: the ME and PH representations are not equal!');
assert(norm(MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))<1e-12, 'FluFluSTD: the ME and PH representations are not equal!');

disp('Test:');
disp('-----');

disp('[alphap, Ap] = FluFluQLD(Qin, Rin, Qout, Rout, true, true):');
[alphap, Ap] = FluFluQLD(Qin, Rin, Qout, Rout, true, true)
disp('[alpha, A] = FluFluQLD(Qin, Rin, Qout, Rout, true):');
[alpha, A] = FluFluQLD(Qin, Rin, Qout, Rout, true)
disp('[betap, Bp] = FluFluSTD(Qin, Rin, Qout, Rout, true, true):');
[betap, Bp] = FluFluSTD(Qin, Rin, Qout, Rout, true, true)
disp('[beta, B] = FluFluSTD(Qin, Rin, Qout, Rout, true):');
[beta, B] = FluFluSTD(Qin, Rin, Qout, Rout, true)

assert(CheckMERepresentation(alpha,A), 'FluFluQLD: invalid ME representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'FluFluSTD: invalid ME representation of the sojourn time!');
assert(CheckPHRepresentation(alphap,Ap), 'FluFluQLD: invalid PH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'FluFluSTD: invalid PH representation of the sojourn time!');

% cross-check
Iin = eye(size(Qin));
Iout = eye(size(Qout));
[gamma, G] = FluidQueueSTD (kron(Qin,Iout)+kron(Iin,Qout), kron(Rin,Iout), kron(Iin,Rout), kron(Qin,Iout)+kron(Rin, pinv(Rout)*Qout));
msmall = MomentsFromME(beta,B,5);
mlarge = MomentsFromME(gamma,G,5);
assert(norm((msmall-mlarge)./msmall)<1e-12, 'FluFluSTD: Large and small model does not give the same results!');

% Little formula
mql = MomentsFromME(alpha,A,1);
mst = MomentsFromME(beta,B,1);
assert(abs(mql-mst*lambda)<1e-12, 'FluFluSTD: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm(MomentsFromPH(alphap,Ap,5)-MomentsFromME(alpha,A,5))<1e-12, 'FluFluQLD: the ME and PH representations are not equal!');
assert(norm(MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))<1e-12, 'FluFluSTD: the ME and PH representations are not equal!');
