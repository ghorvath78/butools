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
help QBDQueue

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

disp('[qld,qlm] = QBDQueue(B, L, F, L0, ''qlDistr'', (0:10), ''qlMoms'', 5):');
[qld,qlm] = QBDQueue(B, L, F, L0, 'qlDistr', (0:10), 'qlMoms', 5)
disp('[std, stm] = QBDQueue(B, L, F, L0, ''stDistr'', (0:0.1:1), ''stMoms'', 5):');
[std, stm] = QBDQueue(B, L, F, L0, 'stDistr', (0:0.1:1), 'stMoms', 5)
disp('[alphap,Ap] = QBDQueue(B, L, F, L0, ''qlDistrDPH''):');
[alphap,Ap] = QBDQueue(B, L, F, L0, 'qlDistrDPH')
disp('[alpha,A] = QBDQueue(B, L, F, L0, ''qlDistrMG''):');
[alpha,A] = QBDQueue(B, L, F, L0, 'qlDistrMG')
disp('[betap,Bp] = QBDQueue(B, L, F, L0, ''stDistrPH''):');
[betap, Bp] = QBDQueue(B, L, F, L0, 'stDistrPH')
disp('[beta,B] = QBDQueue(B, L, F, L0, ''stDistrME''):');
[beta, B] = QBDQueue(B, L, F, L0, 'stDistrME')
 
assert(CheckMGRepresentation(alpha,A), 'QBDQueue: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'QBDQueue: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'QBDQueue: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'QBDQueue: invalid PH representation of the sojourn time!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
assert(abs(mql-mst*lambda)<1e-12, 'QBDQueue: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'QBDQueue: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'QBDQueue: the ME and PH representations are not equal!');

% check moment and distribution calculation
assert(norm(qld-PmfFromMG(alpha,A,(0:10))')<1e-12, 'QBDQueue: qlDistr returns wrong queue length distribution!');
assert(norm(std-CdfFromME(beta,B,(0:0.1:1))')<1e-12, 'QBDQueue: stDistr returns wrong sojourn time distribution!');
assert(norm(qlm-MomentsFromMG(alpha,A,5))<1e-10, 'QBDQueue: qlMoms returns wrong queue length moments!');
assert(norm(stm-MomentsFromME(beta,B,5))<1e-10, 'QBDQueue: stMoms returns wrong sojourn time moments!');

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

disp('[qld,qlm] = QBDQueue(B, L, F, L0, ''qlDistr'', (0:10), ''qlMoms'', 5):');
[qld,qlm] = QBDQueue(B, L, F, L0, 'qlDistr', (0:10), 'qlMoms', 5)
disp('[std, stm] = QBDQueue(B, L, F, L0, ''stDistr'', (0:0.1:1), ''stMoms'', 5):');
[std, stm] = QBDQueue(B, L, F, L0, 'stDistr', (0:0.1:1), 'stMoms', 5)
disp('[alphap,Ap] = QBDQueue(B, L, F, L0, ''qlDistrDPH''):');
[alphap,Ap] = QBDQueue(B, L, F, L0, 'qlDistrDPH')
disp('[alpha,A] = QBDQueue(B, L, F, L0, ''qlDistrMG''):');
[alpha,A] = QBDQueue(B, L, F, L0, 'qlDistrMG')
disp('[betap,Bp] = QBDQueue(B, L, F, L0, ''stDistrPH''):');
[betap, Bp] = QBDQueue(B, L, F, L0, 'stDistrPH')
disp('[beta,B] = QBDQueue(B, L, F, L0, ''stDistrME''):');
[beta, B] = QBDQueue(B, L, F, L0, 'stDistrME')

assert(CheckMGRepresentation(alpha,A), 'QBDQueue: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'QBDQueue: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'QBDQueue: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'QBDQueue: invalid PH representation of the sojourn time!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
assert(abs(mql-mst*lambda)<1e-12, 'QBDQueue: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'QBDQueue: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'QBDQueue: the ME and PH representations are not equal!');

% check moment and distribution calculation
assert(norm(qld-PmfFromMG(alpha,A,(0:10))')<1e-12, 'QBDQueue: qlDistr returns wrong queue length distribution!');
assert(norm(std-CdfFromME(beta,B,(0:0.1:1))')<1e-12, 'QBDQueue: stDistr returns wrong sojourn time distribution!');
assert(norm(qlm-MomentsFromMG(alpha,A,5))<1e-10, 'QBDQueue: qlMoms returns wrong queue length moments!');
assert(norm(stm-MomentsFromME(beta,B,5))<1e-10, 'QBDQueue: stMoms returns wrong sojourn time moments!');

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

disp('[qld,qlm] = QBDQueue(B, L, F, L0, ''qlDistr'', (0:10), ''qlMoms'', 5):');
[qld,qlm] = QBDQueue(B, L, F, L0, 'qlDistr', (0:10), 'qlMoms', 5)
disp('[std, stm] = QBDQueue(B, L, F, L0, ''stDistr'', (0:0.1:1), ''stMoms'', 5):');
[std, stm] = QBDQueue(B, L, F, L0, 'stDistr', (0:0.1:1), 'stMoms', 5)
disp('[alphap,Ap] = QBDQueue(B, L, F, L0, ''qlDistrDPH''):');
[alphap,Ap] = QBDQueue(B, L, F, L0, 'qlDistrDPH')
disp('[alpha,A] = QBDQueue(B, L, F, L0, ''qlDistrMG''):');
[alpha,A] = QBDQueue(B, L, F, L0, 'qlDistrMG')
disp('[betap,Bp] = QBDQueue(B, L, F, L0, ''stDistrPH''):');
[betap, Bp] = QBDQueue(B, L, F, L0, 'stDistrPH')
disp('[beta,B] = QBDQueue(B, L, F, L0, ''stDistrME''):');
[beta, B] = QBDQueue(B, L, F, L0, 'stDistrME')

assert(CheckMGRepresentation(alpha,A), 'QBDQueue: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'QBDQueue: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'QBDQueue: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'QBDQueue: invalid PH representation of the sojourn time!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
assert(abs(mql-mst*lambda)<1e-12, 'QBDQueue: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'QBDQueue: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'QBDQueue: the ME and PH representations are not equal!');

% check moment and distribution calculation
assert(norm(qld-PmfFromMG(alpha,A,(0:10))')<1e-12, 'QBDQueue: qlDistr returns wrong queue length distribution!');
assert(norm(std-CdfFromME(beta,B,(0:0.1:1))')<1e-12, 'QBDQueue: stDistr returns wrong sojourn time distribution!');
assert(norm(qlm-MomentsFromMG(alpha,A,5))<1e-7, 'QBDQueue: qlMoms returns wrong queue length moments!');
assert(norm(stm-MomentsFromME(beta,B,5))<1e-7, 'QBDQueue: stMoms returns wrong sojourn time moments!');

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

disp('[qld,qlm] = QBDQueue(B, L, F, L0, ''qlDistr'', (0:10), ''qlMoms'', 5):');
[qld,qlm] = QBDQueue(B, L, F, L0, 'qlDistr', (0:10), 'qlMoms', 5)
disp('[std, stm] = QBDQueue(B, L, F, L0, ''stDistr'', (0:0.1:1), ''stMoms'', 5):');
[std, stm] = QBDQueue(B, L, F, L0, 'stDistr', (0:0.1:1), 'stMoms', 5)
disp('[alphap,Ap] = QBDQueue(B, L, F, L0, ''qlDistrDPH''):');
[alphap,Ap] = QBDQueue(B, L, F, L0, 'qlDistrDPH')
disp('[alpha,A] = QBDQueue(B, L, F, L0, ''qlDistrMG''):');
[alpha,A] = QBDQueue(B, L, F, L0, 'qlDistrMG')
disp('[betap,Bp] = QBDQueue(B, L, F, L0, ''stDistrPH''):');
[betap, Bp] = QBDQueue(B, L, F, L0, 'stDistrPH')
disp('[beta,B] = QBDQueue(B, L, F, L0, ''stDistrME''):');
[beta, B] = QBDQueue(B, L, F, L0, 'stDistrME')

assert(CheckMGRepresentation(alpha,A), 'QBDQueue: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'QBDQueue: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'QBDQueue: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'QBDQueue: invalid PH representation of the sojourn time!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
assert(abs(mql-mst*lambda)<1e-12, 'QBDQueue: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'QBDQueue: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'QBDQueue: the ME and PH representations are not equal!');

% check moment and distribution calculation
assert(norm(qld-PmfFromMG(alpha,A,(0:10))')<1e-12, 'QBDQueue: qlDistr returns wrong queue length distribution!');
assert(norm(std-CdfFromME(beta,B,(0:0.1:1))')<1e-12, 'QBDQueue: stDistr returns wrong sojourn time distribution!');
assert(norm(qlm-MomentsFromMG(alpha,A,5))<1e-10, 'QBDQueue: qlMoms returns wrong queue length moments!');
assert(norm(stm-MomentsFromME(beta,B,5))<1e-10, 'QBDQueue: stMoms returns wrong sojourn time moments!');

% ============================= MAP/MAP/1 tests ===============================

disp('----------------------------------------------------------------------------');
help MAPMAP1

disp('Input:');
disp('------');

D0=[-8, 2; 1, -3]
D1=[1, 5; 0, 2]

S0=[-10, 4; 0, -7]
S1=[5, 1; 4, 3]

disp('Test:');
disp('-----');

disp('[qld, qlm] = MAPMAP1(D0,D1,S0,S1, ''qlDistr'', (0:10), ''qlMoms'', 5):');
[qld, qlm] = MAPMAP1(D0,D1,S0,S1, 'qlDistr', (0:10), 'qlMoms', 5)
disp('[std, stm] = MAPMAP1(D0,D1,S0,S1, ''stDistr'', (0:0.1:1), ''stMoms'', 5):');
[std, stm] = MAPMAP1(D0,D1,S0,S1, 'stDistr', (0:0.1:1), 'stMoms', 5)
disp('[alphap,Ap] = MAPMAP1(D0,D1,S0,S1,''qlDistrDPH''):');
[alphap,Ap] = MAPMAP1(D0,D1,S0,S1,'qlDistrDPH')
disp('[alpha,A] = MAPMAP1(D0,D1,S0,S1, ''qlDistrMG''):');
[alpha,A] = MAPMAP1(D0,D1,S0,S1, 'qlDistrMG')
disp('[betap,Bp] = MAPMAP1(D0,D1,S0,S1,''stDistrPH''):');
[betap,Bp] = MAPMAP1(D0,D1,S0,S1,'stDistrPH')
disp('[beta,B] = MAPMAP1(D0,D1,S0,S1, ''stDistrME''):');
[beta,B] = MAPMAP1(D0,D1,S0,S1, 'stDistrME')
 
assert(CheckMGRepresentation(alpha,A), 'MAPMAP1: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'MAPMAP1: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'MAPMAP1: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'MAPMAP1: invalid PH representation of the sojourn time!');

% cross-check
IA = eye(size(D0,1));
IS = eye(size(S0,1));
[gamma, G] = QBDQueue (kron(IA,S1), kron(D0,IS)+kron(IA,S0), kron(D1,IS), kron(D0,IS), 'stDistrME');
msmall = MomentsFromME(beta,B,5);
mlarge = MomentsFromME(gamma,G,5);
assert(norm((msmall-mlarge)./msmall)<1e-12, 'MAPMAP1: Large and small model does not give the same results!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
lambda = 1/MarginalMomentsFromMAP(D0,D1,1);
assert(abs(mql-mst*lambda)<1e-12, 'MAPMAP1: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'MAPMAP1: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'MAPMAP1: the ME and PH representations are not equal!');

% check moment and distribution calculation
assert(norm(qld-PmfFromMG(alpha,A,(0:10))')<1e-12, 'MAPMAP1: qlDistr returns wrong queue length distribution!');
assert(norm(std-CdfFromME(beta,B,(0:0.1:1))')<1e-12, 'MAPMAP1: stDistr returns wrong sojourn time distribution!');
assert(norm(qlm-MomentsFromMG(alpha,A,5))<1e-10, 'MAPMAP1: qlMoms returns wrong queue length moments!');
assert(norm(stm-MomentsFromME(beta,B,5))<1e-10, 'MAPMAP1: stMoms returns wrong sojourn time moments!');

disp('Input:');
disp('------');

D0 = [-8, 1, 2; 0, -6, 4; 3 0 -3]
D1 = [4, 1, 0; 0, 2, 0; 0, 0, 0]

disp('Test:');
disp('-----');

disp('[qld, qlm] = MAPMAP1(D0,D1,S0,S1, ''qlDistr'', (0:10), ''qlMoms'', 5):');
[qld, qlm] = MAPMAP1(D0,D1,S0,S1, 'qlDistr', (0:10), 'qlMoms', 5)
disp('[std, stm] = MAPMAP1(D0,D1,S0,S1, ''stDistr'', (0:0.1:1), ''stMoms'', 5):');
[std, stm] = MAPMAP1(D0,D1,S0,S1, 'stDistr', (0:0.1:1), 'stMoms', 5)
disp('[alphap,Ap] = MAPMAP1(D0,D1,S0,S1,''qlDistrDPH''):');
[alphap,Ap] = MAPMAP1(D0,D1,S0,S1,'qlDistrDPH')
disp('[alpha,A] = MAPMAP1(D0,D1,S0,S1, ''qlDistrMG''):');
[alpha,A] = MAPMAP1(D0,D1,S0,S1, 'qlDistrMG')
disp('[betap,Bp] = MAPMAP1(D0,D1,S0,S1,''stDistrPH''):');
[betap,Bp] = MAPMAP1(D0,D1,S0,S1,'stDistrPH')
disp('[beta,B] = MAPMAP1(D0,D1,S0,S1, ''stDistrME''):');
[beta,B] = MAPMAP1(D0,D1,S0,S1, 'stDistrME')
 
assert(CheckMGRepresentation(alpha,A), 'MAPMAP1: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'MAPMAP1: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'MAPMAP1: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'MAPMAP1: invalid PH representation of the sojourn time!');

% cross-check
IA = eye(size(D0,1));
IS = eye(size(S0,1));
[gamma, G] = QBDQueue (kron(IA,S1), kron(D0,IS)+kron(IA,S0), kron(D1,IS), kron(D0,IS), 'stDistrME');
msmall = MomentsFromME(beta,B,5);
mlarge = MomentsFromME(gamma,G,5);
assert(norm((msmall-mlarge)./msmall)<1e-12, 'MAPMAP1: Large and small model does not give the same results!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
lambda = 1/MarginalMomentsFromMAP(D0,D1,1);
assert(abs(mql-mst*lambda)<1e-12, 'MAPMAP1: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'MAPMAP1: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'MAPMAP1: the ME and PH representations are not equal!');

% check moment and distribution calculation
assert(norm(qld-PmfFromMG(alpha,A,(0:10))')<1e-12, 'MAPMAP1: qlDistr returns wrong queue length distribution!');
assert(norm(std-CdfFromME(beta,B,(0:0.1:1))')<1e-12, 'MAPMAP1: stDistr returns wrong sojourn time distribution!');
assert(norm(qlm-MomentsFromMG(alpha,A,5))<1e-10, 'MAPMAP1: qlMoms returns wrong queue length moments!');
assert(norm(stm-MomentsFromME(beta,B,5))<1e-10, 'MAPMAP1: stMoms returns wrong sojourn time moments!');

disp('Input:');
disp('------');

S0 = [-10, 4, 0; 5, -7, 2; 1, 2, -8]
S1 = [0, 0, 6; 0, 0, 0; 0, 3, 2]

D0 = [-8, 1, 2; 0, -6, 4; 3 0 -3]
D1 = [4, 1, 0; 0, 0, 2; 0, 0, 0]

disp('Test:');
disp('-----');

disp('[qld, qlm] = MAPMAP1(D0,D1,S0,S1, ''qlDistr'', (0:10), ''qlMoms'', 5):');
[qld, qlm] = MAPMAP1(D0,D1,S0,S1, 'qlDistr', (0:10), 'qlMoms', 5)
disp('[std, stm] = MAPMAP1(D0,D1,S0,S1, ''stDistr'', (0:0.1:1), ''stMoms'', 5):');
[std, stm] = MAPMAP1(D0,D1,S0,S1, 'stDistr', (0:0.1:1), 'stMoms', 5)
disp('[alphap,Ap] = MAPMAP1(D0,D1,S0,S1,''qlDistrDPH''):');
[alphap,Ap] = MAPMAP1(D0,D1,S0,S1,'qlDistrDPH')
disp('[alpha,A] = MAPMAP1(D0,D1,S0,S1, ''qlDistrMG''):');
[alpha,A] = MAPMAP1(D0,D1,S0,S1, 'qlDistrMG')
disp('[betap,Bp] = MAPMAP1(D0,D1,S0,S1,''stDistrPH''):');
[betap,Bp] = MAPMAP1(D0,D1,S0,S1,'stDistrPH')
disp('[beta,B] = MAPMAP1(D0,D1,S0,S1, ''stDistrME''):');
[beta,B] = MAPMAP1(D0,D1,S0,S1, 'stDistrME')
 
assert(CheckMGRepresentation(alpha,A), 'MAPMAP1: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'MAPMAP1: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'MAPMAP1: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'MAPMAP1: invalid PH representation of the sojourn time!');

% cross-check
IA = eye(size(D0,1));
IS = eye(size(S0,1));
[gamma, G] = QBDQueue (kron(IA,S1), kron(D0,IS)+kron(IA,S0), kron(D1,IS), kron(D0,IS), 'stDistrME');
msmall = MomentsFromME(beta,B,5);
mlarge = MomentsFromME(gamma,G,5);
assert(norm((msmall-mlarge)./msmall)<1e-12, 'MAPMAP1: Large and small model does not give the same results!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
lambda = 1/MarginalMomentsFromMAP(D0,D1,1);
assert(abs(mql-mst*lambda)<1e-12, 'MAPMAP1: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'MAPMAP1: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'MAPMAP1: the ME and PH representations are not equal!');

% check moment and distribution calculation
assert(norm(qld-PmfFromMG(alpha,A,(0:10))')<1e-12, 'MAPMAP1: qlDistr returns wrong queue length distribution!');
assert(norm(std-CdfFromME(beta,B,(0:0.1:1))')<1e-12, 'MAPMAP1: stDistr returns wrong sojourn time distribution!');
assert(norm(qlm-MomentsFromMG(alpha,A,5))<1e-8, 'MAPMAP1: qlMoms returns wrong queue length moments!');
assert(norm(stm-MomentsFromME(beta,B,5))<1e-8, 'MAPMAP1: stMoms returns wrong sojourn time moments!');

% ============================= MAP/PH/1 tests ===============================

disp('----------------------------------------------------------------------------');

disp('Input:');
disp('------');

D0 = [-8, 1, 2; 0, -6, 4; 3 0 -3]
D1 = [4, 1, 0; 0, 0, 2; 0, 0, 0]
sigma = [0.2, 0.7, 0.1]
S = [-10, 4, 0; 5, -7, 2; 1, 2, -8]

S0 = S;
S1 = sum(-S,2)*sigma;

disp('Test:');
disp('-----');

disp('[qld, qlm] = MAPMAP1(D0,D1,S0,S1, ''qlDistr'', (0:10), ''qlMoms'', 5):');
[qld, qlm] = MAPMAP1(D0,D1,S0,S1, 'qlDistr', (0:10), 'qlMoms', 5)
disp('[std, stm] = MAPMAP1(D0,D1,S0,S1, ''stDistr'', (0:0.1:1), ''stMoms'', 5):');
[std, stm] = MAPMAP1(D0,D1,S0,S1, 'stDistr', (0:0.1:1), 'stMoms', 5)
disp('[alphap,Ap] = MAPMAP1(D0,D1,S0,S1,''qlDistrDPH''):');
[alphap,Ap] = MAPMAP1(D0,D1,S0,S1,'qlDistrDPH')
disp('[alpha,A] = MAPMAP1(D0,D1,S0,S1, ''qlDistrMG''):');
[alpha,A] = MAPMAP1(D0,D1,S0,S1, 'qlDistrMG')
disp('[betap,Bp] = MAPMAP1(D0,D1,S0,S1,''stDistrPH''):');
[betap,Bp] = MAPMAP1(D0,D1,S0,S1,'stDistrPH')
disp('[beta,B] = MAPMAP1(D0,D1,S0,S1, ''stDistrME''):');
[beta,B] = MAPMAP1(D0,D1,S0,S1, 'stDistrME')
 
assert(CheckMGRepresentation(alpha,A), 'MAPMAP1: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'MAPMAP1: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'MAPMAP1: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'MAPMAP1: invalid PH representation of the sojourn time!');

% cross-check
IA = eye(size(D0,1));
IS = eye(size(S0,1));
[gamma, G] = QBDQueue (kron(IA,S1), kron(D0,IS)+kron(IA,S0), kron(D1,IS), kron(D0,IS), 'stDistrME');
msmall = MomentsFromME(beta,B,5);
mlarge = MomentsFromME(gamma,G,5);
assert(norm((msmall-mlarge)./msmall)<1e-12, 'MAPMAP1: Large and small model does not give the same results!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
lambda = 1/MarginalMomentsFromMAP(D0,D1,1);
assert(abs(mql-mst*lambda)<1e-12, 'MAPMAP1: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'MAPMAP1: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'MAPMAP1: the ME and PH representations are not equal!');

% check moment and distribution calculation
assert(norm(qld-PmfFromMG(alpha,A,(0:10))')<1e-12, 'MAPMAP1: qlDistr returns wrong queue length distribution!');
assert(norm(std-CdfFromME(beta,B,(0:0.1:1))')<1e-12, 'MAPMAP1: stDistr returns wrong sojourn time distribution!');
assert(norm(qlm-MomentsFromMG(alpha,A,5))<1e-8, 'MAPMAP1: qlMoms returns wrong queue length moments!');
assert(norm(stm-MomentsFromME(beta,B,5))<1e-8, 'MAPMAP1: stMoms returns wrong sojourn time moments!');

% ============================= PH/PH/1 tests ===============================

disp('----------------------------------------------------------------------------');

disp('Input:');
disp('------');

delta = [0.5, 0.1, 0.4]
D = [-8, 1, 2; 0, -6, 4; 3 0 -3]
sigma = [0.2, 0.7, 0.1]
S = [-10, 4, 0; 5, -7, 2; 1, 2, -8]

D0 = D;
D1 = sum(-D,2)*delta;
S0 = S;
S1 = sum(-S,2)*sigma;

disp('Test:');
disp('-----');

disp('[qld, qlm] = MAPMAP1(D0,D1,S0,S1, ''qlDistr'', (0:10), ''qlMoms'', 5):');
[qld, qlm] = MAPMAP1(D0,D1,S0,S1, 'qlDistr', (0:10), 'qlMoms', 5)
disp('[std, stm] = MAPMAP1(D0,D1,S0,S1, ''stDistr'', (0:0.1:1), ''stMoms'', 5):');
[std, stm] = MAPMAP1(D0,D1,S0,S1, 'stDistr', (0:0.1:1), 'stMoms', 5)
disp('[alphap,Ap] = MAPMAP1(D0,D1,S0,S1,''qlDistrDPH''):');
[alphap,Ap] = MAPMAP1(D0,D1,S0,S1,'qlDistrDPH')
disp('[alpha,A] = MAPMAP1(D0,D1,S0,S1, ''qlDistrMG''):');
[alpha,A] = MAPMAP1(D0,D1,S0,S1, 'qlDistrMG')
disp('[betap,Bp] = MAPMAP1(D0,D1,S0,S1,''stDistrPH''):');
[betap,Bp] = MAPMAP1(D0,D1,S0,S1,'stDistrPH')
disp('[beta,B] = MAPMAP1(D0,D1,S0,S1, ''stDistrME''):');
[beta,B] = MAPMAP1(D0,D1,S0,S1, 'stDistrME')
 
assert(CheckMGRepresentation(alpha,A), 'MAPMAP1: invalid MG representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'MAPMAP1: invalid ME representation of the sojourn time!');
assert(CheckDPHRepresentation(alphap,Ap), 'MAPMAP1: invalid DPH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'MAPMAP1: invalid PH representation of the sojourn time!');

% cross-check
IA = eye(size(D0,1));
IS = eye(size(S0,1));
[gamma, G] = QBDQueue (kron(IA,S1), kron(D0,IS)+kron(IA,S0), kron(D1,IS), kron(D0,IS), 'stDistrME');
msmall = MomentsFromME(beta,B,5);
mlarge = MomentsFromME(gamma,G,5);
assert(norm((msmall-mlarge)./msmall)<1e-12, 'MAPMAP1: Large and small model does not give the same results!');

% check Little formula
mql = MomentsFromMG(alpha,A,1);
mst = MomentsFromME(beta,B,1);
lambda = 1/MarginalMomentsFromMAP(D0,D1,1);
assert(abs(mql-mst*lambda)<1e-12, 'MAPMAP1: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm((MomentsFromDPH(alphap,Ap,5)-MomentsFromMG(alpha,A,5))./MomentsFromMG(alpha,A,5))<1e-12, 'MAPMAP1: the MG and DPH representations are not equal!');
assert(norm((MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))./MomentsFromME(beta,B,5))<1e-12, 'MAPMAP1: the ME and PH representations are not equal!');

% check moment and distribution calculation
assert(norm(qld-PmfFromMG(alpha,A,(0:10))')<1e-12, 'MAPMAP1: qlDistr returns wrong queue length distribution!');
assert(norm(std-CdfFromME(beta,B,(0:0.1:1))')<1e-12, 'MAPMAP1: stDistr returns wrong sojourn time distribution!');
assert(norm(qlm-MomentsFromMG(alpha,A,5))<1e-8, 'MAPMAP1: qlMoms returns wrong queue length moments!');
assert(norm(stm-MomentsFromME(beta,B,5))<1e-8, 'MAPMAP1: stMoms returns wrong sojourn time moments!');

% ============================= MMAP[K]/PH[K]/1 PR tests ============================

disp('----------------------------------------------------------------------------');
help MMAPPH1PRPR

disp('Input:');
disp('------');

D0=[-5.49, 0, 1.15, 0; 0, -2.29, 0, 0; 0, 0.08, -1.32, 0; 0.72, 1.17, 0.7, -7.07];
D1=[0.25, 0.38, 0.64, 0; 0, 0, 0, 1.09; 0, 1.24, 0, 0; 0.37, 0, 0, 0];
D2=[0.3, 1.0, 0, 0.48; 0, 0.2, 0, 0; 0, 0, 0, 0; 0.61, 0, 0, 0.2];
D3=[0, 0.98, 0, 0.31; 0, 0, 1.0, 0; 0, 0, 0, 0; 1.1, 0.84, 0.33, 1.03];

sigma3 = [1/6, 5/6];
S3 = [-0.5, 0.5; 0, -3];
sigma2 = [1/1.7, 1-1/1.7];
S2 = [-2/0.85, 2/0.85; 0, -4];
sigma1 = [1/4, 1-1/4];
S1 = [-2.5, 2.5; 0, -10];

disp('Test:');
disp('-----');

[qlm, qld] = MMAPPH1PRPR({D0,D1,D2,D3}, {sigma1,sigma2,sigma3},{S1,S2,S3}, 'qlMoms', 3, 'qlDistr', 500);
momFromDistr = [(0:499); (0:499).^2; (0:499).^3]*qld; 
assert(norm((momFromDistr-qlm)./qlm)<0.001, 'MMAPPH1PRPR: queue length moments and queue length distribution are not consistent!');

[stm, std] = MMAPPH1PRPR({D0,D1,D2,D3}, {sigma1,sigma2,sigma3},{S1,S2,S3}, 'stMoms', 3, 'stDistr', [1,5,10]);

assert(min(min(std))>=0 && max(max(std))<=1 && all(all(diff(std)>=0)), 'MMAPPH1PRPR: invalid sojourn time distribution!');
% check Little formula for every traffic class
lambda1 = 1/MarginalMomentsFromMAP(D0+D2+D3,D1,1);
lambda2 = 1/MarginalMomentsFromMAP(D0+D1+D3,D2,1);
lambda3 = 1/MarginalMomentsFromMAP(D0+D1+D2,D3,1);
assert(abs(qlm(1,1)-stm(1,1)*lambda1)<1e-12, 'MMAPPH1PRPR: Little formula does not hold for class 1!');
assert(abs(qlm(1,2)-stm(1,2)*lambda2)<1e-12, 'MMAPPH1PRPR: Little formula does not hold for class 2!');
assert(abs(qlm(1,3)-stm(1,3)*lambda3)<1e-12, 'MMAPPH1PRPR: Little formula does not hold for class 3!');

% ============================= MMAP[K]/PH[K]/1 NP tests ============================

disp('----------------------------------------------------------------------------');
help MMAPPH1NPPR

disp('Input:');
disp('------');

D0=[-5.49, 0, 1.15, 0; 0, -2.29, 0, 0; 0, 0.08, -1.32, 0; 0.72, 1.17, 0.7, -7.07];
D1=[0.25, 0.38, 0.64, 0; 0, 0, 0, 1.09; 0, 1.24, 0, 0; 0.37, 0, 0, 0];
D2=[0.3, 1.0, 0, 0.48; 0, 0.2, 0, 0; 0, 0, 0, 0; 0.61, 0, 0, 0.2];
D3=[0, 0.98, 0, 0.31; 0, 0, 1.0, 0; 0, 0, 0, 0; 1.1, 0.84, 0.33, 1.03];

sigma3 = [1/6, 5/6];
S3 = [-0.5, 0.5; 0, -3];
sigma2 = [1/1.7, 1-1/1.7];
S2 = [-2/0.85, 2/0.85; 0, -4];
sigma1 = [1/4, 1-1/4];
S1 = [-2.5, 2.5; 0, -10];

disp('Test:');
disp('-----');

[qlm, qld] = MMAPPH1NPPR({D0,D1,D2,D3}, {sigma1,sigma2,sigma3},{S1,S2,S3}, 'qlMoms', 3, 'qlDistr', 500);
momFromDistr = [(0:499); (0:499).^2; (0:499).^3]*qld; 
assert(norm((momFromDistr-qlm)./qlm)<0.001, 'MMAPPH1PRPR: queue length moments and queue length distribution are not consistent!');

[stm, std] = MMAPPH1NPPR({D0,D1,D2,D3}, {sigma1,sigma2,sigma3},{S1,S2,S3}, 'stMoms', 3, 'stDistr', [1,5,10]);

assert(min(min(std))>=0 && max(max(std))<=1 && all(all(diff(std)>=0)), 'MMAPPH1NPPR: invalid sojourn time distribution!');
% check Little formula for every traffic class
lambda1 = 1/MarginalMomentsFromMAP(D0+D2+D3,D1,1);
lambda2 = 1/MarginalMomentsFromMAP(D0+D1+D3,D2,1);
lambda3 = 1/MarginalMomentsFromMAP(D0+D1+D2,D3,1);
assert(abs(qlm(1,1)-stm(1,1)*lambda1)<1e-12, 'MMAPPH1NPPR: Little formula does not hold for class 1!');
assert(abs(qlm(1,2)-stm(1,2)*lambda2)<1e-12, 'MMAPPH1NPPR: Little formula does not hold for class 2!');
assert(abs(qlm(1,3)-stm(1,3)*lambda3)<1e-12, 'MMAPPH1NPPR: Little formula does not hold for class 3!');


% ============================= Fluid tests ===============================

disp('----------------------------------------------------------------------------');
help FluidQueue

disp('Input:');
disp('------');
Q = [-9 2 4 0 1 2; 6 -25 5 3 7 4; 1 3 -4 0 0 0; 0 0 0 -8 3 5; 7 3 0 2 -13 1; 7 8 0 3 8 -26]
Rin = diag([4 2 1 0 0 3])
Rout = diag([6 2 0 0 3 2])
lambda = sum(CTMCSolve(Q)*Rin);

disp('Test:');
disp('-----');

disp('[qld, qlm] = FluidQueue(Q, Rin, Rout, ''qlDistr'', (0:0.1:1), ''qlMoms'', 5):');
[qld, qlm] = FluidQueue(Q, Rin, Rout, 'qlDistr', (0:0.1:1), 'qlMoms', 5)
disp('[std, stm] = FluidQueue(Q, Rin, Rout, ''stDistr'', (0:0.1:1), ''stMoms'', 5):');
[std, stm] = FluidQueue(Q, Rin, Rout, 'stDistr', (0:0.1:1), 'stMoms', 5)
disp('[alpha, A] = FluidQueue(Q, Rin, Rout, ''qlDistrME''):');
[alpha, A] = FluidQueue(Q, Rin, Rout, 'qlDistrME')
disp('[alphap, Ap] = FluidQueue(Q, Rin, Rout, ''qlDistrPH''):');
[alphap, Ap] = FluidQueue(Q, Rin, Rout, 'qlDistrPH')
disp('[beta, B] = FluidQueue(Q, Rin, Rout, ''stDistrME''):');
[beta, B] = FluidQueue(Q, Rin, Rout, 'stDistrME')
disp('[betap, Bp] = FluidQueue(Q, Rin, Rout, ''stDistrPH''):');
[betap, Bp] = FluidQueue(Q, Rin, Rout, 'stDistrPH')

assert(CheckMERepresentation(alpha,A), 'FluidQueue: invalid ME representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'FluidQueue: invalid ME representation of the sojourn time!');
assert(CheckPHRepresentation(alphap,Ap), 'FluidQueue: invalid PH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'FluidQueue: invalid PH representation of the sojourn time!');

% check Little formula
mql = MomentsFromME(alpha,A,1);
mst = MomentsFromME(beta,B,1);
assert(abs(mql-mst*lambda)<1e-12, 'FluidQueue: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm(MomentsFromPH(alphap,Ap,5)-MomentsFromME(alpha,A,5))<1e-12, 'FluidQueue: the ME and PH representations are not equal!');
assert(norm(MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))<1e-12, 'FluidQueue: the ME and PH representations are not equal!');

% check moment and distribution calculation
assert(norm(qld-CdfFromME(alpha,A,(0:0.1:1))')<1e-12, 'FluidQueue: qlDistr returns wrong fluid level distribution!');
assert(norm(std-CdfFromME(beta,B,(0:0.1:1))')<1e-12, 'FluidQueue: stDistr returns wrong sojourn time distribution!');
assert(norm(qlm-MomentsFromME(alpha,A,5))<1e-8, 'FluidQueue: qlMoms returns wrong fluid level moments!');
assert(norm(stm-MomentsFromME(beta,B,5))<1e-8, 'FluidQueue: stMoms returns wrong sojourn time moments!');

% ============================= FluFlu tests ===============================

disp('----------------------------------------------------------------------------');
help FluFluQueue

disp('Input:');
disp('------');

Qin = [-2 1 1; 2 -5 3; 4 0 -4]
Rin = diag([3 7 0])
lambda = sum(CTMCSolve(Qin)*Rin);

Qout = [-4 1 3; 6 -8 2; 3 7 -10]
Rout = diag([1 7 15])

disp('Test:');
disp('-----');

disp('[qld, qlm] = FluFluQueue(Qin, Rin, Qout, Rout, false, ''qlDistr'', (0:0.1:1), ''qlMoms'', 5):');
[qld, qlm] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'qlDistr', (0:0.1:1), 'qlMoms', 5)
disp('[std, stm] = FluFluQueue(Qin, Rin, Qout, Rout, false, ''stDistr'', (0:0.1:1), ''stMoms'', 5):');
[std, stm] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'stDistr', (0:0.1:1), 'stMoms', 5)
disp('[alpha, A] = FluFluQueue(Qin, Rin, Qout, Rout, false, ''qlDistrME''):');
[alpha, A] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'qlDistrME')
disp('[alphap, Ap] = FluFluQueue(Qin, Rin, Qout, Rout, false, ''qlDistrPH''):');
[alphap, Ap] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'qlDistrPH')
disp('[beta, B] = FluFluQueue(Qin, Rin, Qout, Rout, false, ''stDistrME''):');
[beta, B] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'stDistrME')
disp('[betap, Bp] = FluFluQueue(Qin, Rin, Qout, Rout, false, ''stDistrPH''):');
[betap, Bp] = FluFluQueue(Qin, Rin, Qout, Rout, false, 'stDistrPH')

assert(CheckMERepresentation(alpha,A), 'FluFluQueue: invalid ME representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'FluFluQueue: invalid ME representation of the sojourn time!');
assert(CheckPHRepresentation(alphap,Ap), 'FluFluQueue: invalid PH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'FluFluQueue: invalid PH representation of the sojourn time!');

% cross-check
Iin = eye(size(Qin));
Iout = eye(size(Qout));
[gamma, G] = FluidQueue (kron(Qin,Iout)+kron(Iin,Qout), kron(Rin,Iout), kron(Iin,Rout), 'stDistrME');
msmall = MomentsFromME(beta,B,5);
mlarge = MomentsFromME(gamma,G,5);
assert(norm((msmall-mlarge)./msmall)<1e-12, 'FluFluQueue: Large and small model does not give the same results!');

% Little formula
mql = MomentsFromME(alpha,A,1);
mst = MomentsFromME(beta,B,1);
assert(abs(mql-mst*lambda)<1e-12, 'FluFluQueue: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm(MomentsFromPH(alphap,Ap,5)-MomentsFromME(alpha,A,5))<1e-12, 'FluFluQueue: the ME and PH representations are not equal!');
assert(norm(MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))<1e-12, 'FluFluQueue: the ME and PH representations are not equal!');

% check moment and distribution calculation
assert(norm(qld-CdfFromME(alpha,A,(0:0.1:1))')<1e-12, 'FluFluQueue: qlDistr returns wrong fluid level distribution!');
assert(norm(std-CdfFromME(beta,B,(0:0.1:1))')<1e-12, 'FluFluQueue: stDistr returns wrong sojourn time distribution!');
assert(norm(qlm-MomentsFromME(alpha,A,5))<1e-8, 'FluFluQueue: qlMoms returns wrong fluid level moments!');
assert(norm(stm-MomentsFromME(beta,B,5))<1e-8, 'FluFluQueue: stMoms returns wrong sojourn time moments!');

disp('Test:');
disp('-----');

disp('[qld, qlm] = FluFluQueue(Qin, Rin, Qout, Rout, true, ''qlDistr'', (0:0.1:1), ''qlMoms'', 5):');
[qld, qlm] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'qlDistr', (0:0.1:1), 'qlMoms', 5)
disp('[std, stm] = FluFluQueue(Qin, Rin, Qout, Rout, true, ''stDistr'', (0:0.1:1), ''stMoms'', 5):');
[std, stm] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'stDistr', (0:0.1:1), 'stMoms', 5)
disp('[alpha, A] = FluFluQueue(Qin, Rin, Qout, Rout, true, ''qlDistrME''):');
[alpha, A] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'qlDistrME')
disp('[alphap, Ap] = FluFluQueue(Qin, Rin, Qout, Rout, true, ''qlDistrPH''):');
[alphap, Ap] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'qlDistrPH')
disp('[beta, B] = FluFluQueue(Qin, Rin, Qout, Rout, true, ''stDistrME''):');
[beta, B] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'stDistrME')
disp('[betap, Bp] = FluFluQueue(Qin, Rin, Qout, Rout, true, ''stDistrPH''):');
[betap, Bp] = FluFluQueue(Qin, Rin, Qout, Rout, true, 'stDistrPH')

assert(CheckMERepresentation(alpha,A), 'FluFluQueue: invalid ME representation of the queue length!');
assert(CheckMERepresentation(beta,B), 'FluFluQueue: invalid ME representation of the sojourn time!');
assert(CheckPHRepresentation(alphap,Ap), 'FluFluQueue: invalid PH representation of the queue length!');
assert(CheckPHRepresentation(betap,Bp), 'FluFluQueue: invalid PH representation of the sojourn time!');

% cross-check
Iin = eye(size(Qin));
Iout = eye(size(Qout));
[gamma, G] = FluidQueue (kron(Qin,Iout)+kron(Iin,Qout), kron(Rin,Iout), kron(Iin,Rout), 'Q0', kron(Qin,Iout)+kron(Rin,pinv(Rout)*Qout), 'stDistrME');
msmall = MomentsFromME(beta,B,5);
mlarge = MomentsFromME(gamma,G,5);
assert(norm((msmall-mlarge)./msmall)<1e-12, 'FluFluQueue: Large and small model does not give the same results!');

% Little formula
mql = MomentsFromME(alpha,A,1);
mst = MomentsFromME(beta,B,1);
assert(abs(mql-mst*lambda)<1e-12, 'FluFluQueue: Little formula does not hold!');

% check the equality of the PH and ME results
assert(norm(MomentsFromPH(alphap,Ap,5)-MomentsFromME(alpha,A,5))<1e-12, 'FluFluQueue: the ME and PH representations are not equal!');
assert(norm(MomentsFromPH(betap,Bp,5)-MomentsFromME(beta,B,5))<1e-12, 'FluFluQueue: the ME and PH representations are not equal!');

% check moment and distribution calculation
assert(norm(qld-CdfFromME(alpha,A,(0:0.1:1))')<1e-12, 'FluFluQueue: qlDistr returns wrong fluid level distribution!');
assert(norm(std-CdfFromME(beta,B,(0:0.1:1))')<1e-12, 'FluFluQueue: stDistr returns wrong sojourn time distribution!');
assert(norm(qlm-MomentsFromME(alpha,A,5))<1e-8, 'FluFluQueue: qlMoms returns wrong fluid level moments!');
assert(norm(stm-MomentsFromME(beta,B,5))<1e-8, 'FluFluQueue: stMoms returns wrong sojourn time moments!');
