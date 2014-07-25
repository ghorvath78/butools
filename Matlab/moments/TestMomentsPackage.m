format short g;

disp('---BuTools: Moments package test file---');

disp('Enable the verbose messages with the BuToolsVerbose flag')
global BuToolsVerbose;
BuToolsVerbose = true

disp('Enable input parameter checking with the BuToolsCheckInput flag')
global BuToolsCheckInput;
BuToolsCheckInput = true

disp('----------------------------------------------------------------------------');
help NormMomsFromMoms
help MomsFromNormMoms

disp('Test:');
disp('-----');

M = [1.2 5 38 495 9215]
disp('nmoms=NormMomsFromMoms(M)');
nmoms=NormMomsFromMoms(M);
disp(nmoms);
disp('moms=MomsFromNormMoms(nmoms)');
moms=MomsFromNormMoms(nmoms);
disp(moms);
assert(max(abs(moms-M))<1e-14, 'Calling the moment conversion and its inverse did not give back the original moments!');

disp('----------------------------------------------------------------------------');
help ReducedMomsFromMoms
help MomsFromReducedMoms

disp('Test:');
disp('-----');

M
disp('rmoms=ReducedMomsFromMoms(M)');
rmoms=ReducedMomsFromMoms(M);
disp(rmoms);
disp('moms=MomsFromReducedMoms(rmoms)');
moms=MomsFromReducedMoms(rmoms);
disp(moms);
assert(max(abs(moms-M))<1e-14, 'Calling the moment conversion and its inverse did not give back the original moments!');

disp('----------------------------------------------------------------------------');
help FactorialMomsFromMoms
help MomsFromFactorialMoms

disp('Test:');
disp('-----');

M = [1.3 2.4 6.03 20.5 89.5 474.9]
disp('fmoms=FactorialMomsFromMoms(M)');
fmoms=FactorialMomsFromMoms(M);
disp(fmoms);
disp('moms=MomsFromFactorialmoms(fmoms)');
moms=MomsFromFactorialMoms(fmoms);
disp(moms)
assert(max(abs(moms-M))<1e-14, 'Calling the moment conversion and its inverse did not give back the original moments!');


disp('----------------------------------------------------------------------------');
help JFactorialMomsFromJMoms
help JMomsFromJFactorialMoms

disp('Test:');
disp('-----');

N = [0.7 2 3 4; 5 6 7 8; 9 10 11 12]
disp('JFmoms=JFactorialMomsFromJMoms(N)');
JFmoms=JFactorialMomsFromJMoms(N);
disp(JFmoms);
disp('Jmoms=JMomsFromJFactorialMoms(JFmoms)');
Jmoms=JMomsFromJFactorialMoms(JFmoms);
disp(Jmoms);
assert(max(max(abs(moms-M)))<1e-14, 'Calling the moment conversion and its inverse did not give back the original moments!');

disp('----------------------------------------------------------------------------');
help HankelMomsFromMoms
help MomsFromHankelMoms

disp('Test:');
disp('-----');

M
disp('hmoms=HankelMomsFromMoms(M)');
hmoms=HankelMomsFromMoms(M);
disp(hmoms);
disp('moms=MomsFromHankelMoms(hmoms)');
moms=MomsFromHankelMoms(hmoms);
disp(moms);
assert(max(abs(moms-M))<1e-12, 'Calling the moment conversion and its inverse did not give back the original moments!');

disp('----------------------------------------------------------------------------');

help CheckMoments

disp('Test:');
disp('-----');

M = [1.2 5 8 29 3412]
disp('flag=CheckMoments(M)');
flag=CheckMoments(M);
disp(flag);
assert(flag==0, 'CheckMoments did not recognize a valid moment sequence!');


M = [1.3 2.4 6.03 20.5 89.5]
disp('flag=CheckMoments(M)');
flag=CheckMoments(M);
disp(flag);
assert(flag==1, 'CheckMoments did not recognize an invalid moment sequence!');
