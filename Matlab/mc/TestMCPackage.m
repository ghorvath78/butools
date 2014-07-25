format short g;

disp('---BuTools: MC package test file---');

disp('Enable the verbose messages with the BuToolsVerbose flag')
global BuToolsVerbose;
BuToolsVerbose = true

disp('Enable input parameter checking with the BuToolsCheckInput flag')
global BuToolsCheckInput;
BuToolsCheckInput = true

disp('----------------------------------------------------------------------------');
help CRPSolve

disp('Input:');
disp('------');

Q1= [ -4.3 3.5 0.8; -8.4 6.5 1.9; 17.3 -12.7 -4.6]

disp('Test:');
disp('-----');

disp('CRPSolve(Q1):');
ret=CRPSolve(Q1);
disp(ret);
assert(abs(max(ret*Q1))<1e-14, 'The solution does not satisfy ret*Q1=0!');

disp('----------------------------------------------------------------------------');
help DRPSolve

disp('Input:');
disp('------');

Q2=[-0.9 0.5 1.4; 0.9 -0.9 1; 0.3 1.3 -0.6]

disp('Test:');
disp('-----');

disp('DRPSolve(Q2):');
ret=DRPSolve(Q2);
disp(ret);
assert(abs(max(ret*Q2-ret))<1e-14, 'The solution does not satisfy ret*Q2=ret!');

disp('----------------------------------------------------------------------------');
help CTMCSolve

disp('Input:');
disp('------');

Q3 = [0.1 0.5 0.4; 0.9 0.1 0]

disp('Test:');
disp('-----');

disp('CTMCSolve(Q3):');
try
    ret=CTMCSolve(Q3);
    disp(ret);
catch err
    disp(err);
end

disp('Input:');
disp('------');

Q4 = [0.1 0.5 0.4; 0.9 0.1 0; 0.9 0.1 0]

disp('Test:');
disp('-----');

disp('CTMCSolve(Q4):');

try
    ret=CTMCSolve(q);
    disp(ret);
catch err
    disp(err);
end
    
disp('Input:');
disp('------');

Q5 = [-0.9 0.5 0.4; 0.9 -0.9 0; 0.3 0.3 -0.6]

disp('Test:');
disp('-----');

disp('CTMCSolve(Q5):');
ret=CTMCSolve(Q5);
disp(ret);
assert(abs(max(ret*Q5))<1e-14, 'The solution does not satisfy ret*Q5=0!');

disp('----------------------------------------------------------------------------');
help DTMCSolve

disp('Input:');
disp('------');

Q6 = [0.1 0.5 0.4; 0.9 0.1 0; 0.3 -0.3 0.4]

disp('Test:');
disp('-----');

disp('DTMCSolve(Q6):');
try
    ret=DTMCSolve(Q6);
    disp(ret);
catch err
    disp(err);
end
    
disp('Input:');
disp('------');

Q7 = [0.1 0.5 0.4; 0.9 0.1 0; 0.3 0.3 0.4]

disp('Test:');
disp('-----');

disp('DTMCSolve(Q7):');
ret=DTMCSolve(Q7);
disp(ret);
assert(abs(max(ret*Q7-ret))<1e-14, 'The solution does not satisfy ret*Q7=ret!');

disp('----------------------------------------------------------------------------');
help CheckGenerator

disp('Input:');
disp('------');

Q8 = [-0.9 0.2 0.4; 0 0.9 0.9; 0 0.6 -0.6]

disp('Test:');
disp('-----');

disp('CheckGenerator(Q8,true):');
flag=CheckGenerator(Q8,true);
disp(flag);
assert(flag==0,'CheckGenerator did not detect bad row sum!');

disp('Input:');
disp('------');

Q9 = [-0.9 0.5 0.4; 0.9 -0.9 0; 0.3 0.3 -0.6]

disp('Test:');
disp('-----');

disp('CheckGenerator(Q9,true):');
flag=CheckGenerator(Q9,true);
disp(flag);
assert(flag==1,'CheckGenerator did not recognize a valid input!');


disp('Input:');
disp('------');

Q10 = [-0.9 0.2 0.4; 0.9 -0.9 0; 0.3 0.3 -0.6]

disp('Test:');
disp('-----');

disp('CheckGenerator(Q10,true):');
flag=CheckGenerator(Q10,true);
disp(flag);
assert(flag==1,'CheckGenerator did not recognize a valid input!');

disp('Input:');
disp('------');

Q11 = [-0.9 0.5 0.4; 0.9 -1.1 0; 0.3 0.3 -0.6]

disp('Test:');
disp('-----');

disp('CheckGenerator(Q11):');
flag=CheckGenerator(Q11);
disp(flag);
assert(flag==0,'CheckGenerator did not recognize the non-zero row sum!');

disp('Input:');
disp('------');

Q12 = [-0.9 0.5 0.4; 0.9 -0.9 0; 0.3 0.3 -0.6]

disp('Test:');
disp('-----');

disp('CheckGenerator(Q12):');
flag=CheckGenerator(Q12);
disp(flag);
assert(flag==1,'CheckGenerator did not recognize a valid input!');

disp('-------------------------------------------------------------------------');
help CheckProbMatrix

disp('Input:');
disp('------');

Q13 = [0.1 0.5 0.4; 0.9 0.1 0; 0.3 -0.1 0.4]

disp('Test:');
disp('-----');

disp('CheckProbMatrix(Q13):');
flag=CheckProbMatrix(Q13);
disp(flag);
assert(flag==0,'CheckProbMatrix did not recognize the negative entry!');

disp('Input:');
disp('------');

Q14 = [0.1 0.5 0.4; 0.9 0.1 0; 0.3 0.1 0.4]

disp('Test:');
disp('-----');

disp('CheckProbMatrix(Q14):');
flag=CheckProbMatrix(Q14);
disp(flag);
assert(flag==0,'CheckProbMatrix did not recognize the invalid row sum!');

disp('Input:');
disp('------');

Q15 = [0.1 0.5 0.4; 0.9 0.1 0; 0.3 0.3 0.4]

disp('Test:');
disp('-----');

disp('CheckProbMatrix(Q15):');
flag=CheckProbMatrix(Q15);
disp(flag);
assert(flag==1,'CheckProbMatrix did not recognize that the input is valid!');

disp('Input:');
disp('------');

Q16 = [0.1 0.5 0.4; 0.9 0.1 0; 0.3 0.3 0.4]

disp('Test:');
disp('-----');

disp('CheckProbMatrix(Q16,true):');
flag=CheckProbMatrix(Q16,true);
disp(flag);
assert(flag==0,'CheckProbMatrix did not recognize wrong transient matrix!');

disp('Input:');
disp('------');

Q17 = [0.1 0.5 0.4; 0.9 0.1 0; 0.3 0.1 0.4]

disp('Test:');
disp('-----');

disp('CheckProbMatrix(Q17,true):');
flag=CheckProbMatrix(Q17,true);
disp(flag);
assert(flag==1,'CheckProbMatrix did not recognize that the input is valid!');

disp('-------------------------------------------------------------------------');
help CheckProbVector

disp('Input:');
disp('------');

Q18 = [1.1 -0.1]

disp('Test:');
disp('-----');

disp('CheckProbVector(Q18):');
flag=CheckProbVector(Q18);
disp(flag);
assert(flag==0,'CheckProbVector did not recognize the negative entry!');

disp('Input:');
disp('------');

Q19 = [1.1 0.1]

disp('Test:');
disp('-----');

disp('CheckProbVector(Q19):');
flag=CheckProbVector(Q19);
disp(flag);
assert(flag==0,'CheckProbVector did not recognize invalid sum!');

disp('Input:');
disp('------');

Q20 = [1 0]

disp('Test:');
disp('-----');

disp('CheckProbVector(Q20):');
flag=CheckProbVector(Q20);
disp(flag);
assert(flag==1,'CheckProbVector did not recognize that the input is valid!');

disp('Input:');
disp('------');

Q21 = [0.9 -0.1]

disp('Test:');
disp('-----');

disp('CheckProbVector(Q21,true):');
flag=CheckProbVector(Q21,true);
disp(flag);
assert(flag==0,'CheckProbVector did not recognize the negative entry!');

disp('Input:');
disp('------');

Q22 = [0.9 0.1]

disp('Test:');
disp('-----');

disp('CheckProbVector(Q22,true):');
flag=CheckProbVector(Q22,true);
disp(flag);
assert(flag==1,'CheckProbVector did not recognize that the prob. vector is not transient!');

disp('Input:');
disp('------');

Q23 = [0.8 0.1]

disp('Test:');
disp('-----');

disp('CheckProbVector(Q23,true):');
flag=CheckProbVector(Q23,true);
disp(flag);
assert(flag==1,'CheckProbVector did not recognize that the input is valid!');
