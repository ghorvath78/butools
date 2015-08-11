clear all
run('/home/gabor/github/butools/Matlab/BuToolsInit.m')

disp('---BuTools: Moments package test file---');
disp('Enable the verbose messages with the BuToolsVerbose flag');
global BuToolsVerbose;
BuToolsVerbose = true;
disp('Enable input parameter checking with the BuToolsCheckInput flag');
global BuToolsCheckInput;
BuToolsCheckInput = true;
global BuToolsCheckPrecision;
format short g;
disp('========================================')
disp('Testing BuTools function NormMomsFromMoms')
disp('Input:');
disp('------');
M = [1.2, 5., 38., 495., 9215.];
disp('M = ');
disp(M);
disp('Test:');
disp('-----');
disp('nmoms = NormMomsFromMoms(M);:');
nmoms = NormMomsFromMoms(M);
disp('nmoms = ');
disp(nmoms);
disp('moms = MomsFromNormMoms(nmoms);:');
moms = MomsFromNormMoms(nmoms);
disp('moms = ');
disp(moms);
disp('err = norm(moms-M);:');
err = norm(moms-M);
disp('err = ');
disp(err);
assert(err<10^-10, 'Calling the moment conversion and its inverse did not give back the original moments!');
disp('========================================')
disp('Testing BuTools function MomsFromNormMoms')
disp('Input:');
disp('------');
M = [1.2, 5., 38., 495., 9215.];
disp('M = ');
disp(M);
disp('Test:');
disp('-----');
disp('nmoms = NormMomsFromMoms(M);:');
nmoms = NormMomsFromMoms(M);
disp('nmoms = ');
disp(nmoms);
disp('moms = MomsFromNormMoms(nmoms);:');
moms = MomsFromNormMoms(nmoms);
disp('moms = ');
disp(moms);
disp('err = norm(moms-M);:');
err = norm(moms-M);
disp('err = ');
disp(err);
assert(err<10^-10, 'Calling the moment conversion and its inverse did not give back the original moments!');
disp('========================================')
disp('Testing BuTools function ReducedMomsFromMoms')
disp('Input:');
disp('------');
M = [1.2, 5., 38., 495., 9215.];
disp('M = ');
disp(M);
disp('Test:');
disp('-----');
disp('rmoms = ReducedMomsFromMoms(M);:');
rmoms = ReducedMomsFromMoms(M);
disp('rmoms = ');
disp(rmoms);
disp('moms = MomsFromReducedMoms(rmoms);:');
moms = MomsFromReducedMoms(rmoms);
disp('moms = ');
disp(moms);
disp('err = norm(moms-M);:');
err = norm(moms-M);
disp('err = ');
disp(err);
assert(err<10^-10, 'Calling the moment conversion and its inverse did not give back the original moments!');
disp('========================================')
disp('Testing BuTools function MomsFromReducedMoms')
disp('Input:');
disp('------');
M = [1.2, 5., 38., 495., 9215.];
disp('M = ');
disp(M);
disp('Test:');
disp('-----');
disp('rmoms = ReducedMomsFromMoms(M);:');
rmoms = ReducedMomsFromMoms(M);
disp('rmoms = ');
disp(rmoms);
disp('moms = MomsFromReducedMoms(rmoms);:');
moms = MomsFromReducedMoms(rmoms);
disp('moms = ');
disp(moms);
disp('err = norm(moms-M);:');
err = norm(moms-M);
disp('err = ');
disp(err);
assert(err<10^-10, 'Calling the moment conversion and its inverse did not give back the original moments!');
disp('========================================')
disp('Testing BuTools function FactorialMomsFromMoms')
disp('Input:');
disp('------');
M = [1.3, 2.4, 6.03, 20.5, 89.5, 474.9];
disp('M = ');
disp(M);
disp('Test:');
disp('-----');
disp('fmoms = FactorialMomsFromMoms(M);:');
fmoms = FactorialMomsFromMoms(M);
disp('fmoms = ');
disp(fmoms);
disp('moms = MomsFromFactorialMoms(fmoms);:');
moms = MomsFromFactorialMoms(fmoms);
disp('moms = ');
disp(moms);
disp('err = norm(moms-M);:');
err = norm(moms-M);
disp('err = ');
disp(err);
assert(err<10^-10, 'Calling the moment conversion and its inverse did not give back the original moments!');
disp('========================================')
disp('Testing BuTools function MomsFromFactorialMoms')
disp('Input:');
disp('------');
M = [1.3, 2.4, 6.03, 20.5, 89.5, 474.9];
disp('M = ');
disp(M);
disp('Test:');
disp('-----');
disp('fmoms = FactorialMomsFromMoms(M);:');
fmoms = FactorialMomsFromMoms(M);
disp('fmoms = ');
disp(fmoms);
disp('moms = MomsFromFactorialMoms(fmoms);:');
moms = MomsFromFactorialMoms(fmoms);
disp('moms = ');
disp(moms);
disp('err = norm(moms-M);:');
err = norm(moms-M);
disp('err = ');
disp(err);
assert(err<10^-10, 'Calling the moment conversion and its inverse did not give back the original moments!');
disp('========================================')
disp('Testing BuTools function JFactorialMomsFromJMoms')
disp('Input:');
disp('------');
MM = [0.7, 2., 3., 4.; 5., 6., 7., 8.; 9., 10., 11., 12.];
disp('MM = ');
disp(MM);
disp('Test:');
disp('-----');
disp('JFmoms = JFactorialMomsFromJMoms(MM);:');
JFmoms = JFactorialMomsFromJMoms(MM);
disp('JFmoms = ');
disp(JFmoms);
disp('Jmoms = JMomsFromJFactorialMoms(JFmoms);:');
Jmoms = JMomsFromJFactorialMoms(JFmoms);
disp('Jmoms = ');
disp(Jmoms);
disp('err = norm(Jmoms-MM);:');
err = norm(Jmoms-MM);
disp('err = ');
disp(err);
assert(err<10^-10, 'Calling the moment conversion and its inverse did not give back the original moments!');
disp('========================================')
disp('Testing BuTools function JMomsFromJFactorialMoms')
disp('Input:');
disp('------');
MM = [0.7, 2., 3., 4.; 5., 6., 7., 8.; 9., 10., 11., 12.];
disp('MM = ');
disp(MM);
disp('Test:');
disp('-----');
disp('JFmoms = JFactorialMomsFromJMoms(MM);:');
JFmoms = JFactorialMomsFromJMoms(MM);
disp('JFmoms = ');
disp(JFmoms);
disp('Jmoms = JMomsFromJFactorialMoms(JFmoms);:');
Jmoms = JMomsFromJFactorialMoms(JFmoms);
disp('Jmoms = ');
disp(Jmoms);
disp('err = norm(Jmoms-MM);:');
err = norm(Jmoms-MM);
disp('err = ');
disp(err);
assert(err<10^-10, 'Calling the moment conversion and its inverse did not give back the original moments!');
disp('========================================')
disp('Testing BuTools function HankelMomsFromMoms')
disp('Input:');
disp('------');
M = [1.3, 2.4, 6.03, 20.5, 89.5, 474.9];
disp('M = ');
disp(M);
disp('Test:');
disp('-----');
disp('hmoms = HankelMomsFromMoms(M);:');
hmoms = HankelMomsFromMoms(M);
disp('hmoms = ');
disp(hmoms);
disp('moms = MomsFromHankelMoms(hmoms);:');
moms = MomsFromHankelMoms(hmoms);
disp('moms = ');
disp(moms);
disp('err = norm(moms-M);:');
err = norm(moms-M);
disp('err = ');
disp(err);
assert(err<10^-10, 'Calling the moment conversion and its inverse did not give back the original moments!');
disp('========================================')
disp('Testing BuTools function MomsFromHankelMoms')
disp('Input:');
disp('------');
M = [1.3, 2.4, 6.03, 20.5, 89.5, 474.9];
disp('M = ');
disp(M);
disp('Test:');
disp('-----');
disp('hmoms = HankelMomsFromMoms(M);:');
hmoms = HankelMomsFromMoms(M);
disp('hmoms = ');
disp(hmoms);
disp('moms = MomsFromHankelMoms(hmoms);:');
moms = MomsFromHankelMoms(hmoms);
disp('moms = ');
disp(moms);
disp('err = norm(moms-M);:');
err = norm(moms-M);
disp('err = ');
disp(err);
assert(err<10^-10, 'Calling the moment conversion and its inverse did not give back the original moments!');
disp('========================================')
disp('Testing BuTools function CheckMoments')
disp('Input:');
disp('------');
M = [1.2, 5., 8., 29., 3412.];
disp('M = ');
disp(M);
disp('Test:');
disp('-----');
disp('flag = CheckMoments(M);:');
flag = CheckMoments(M);
disp('flag = ');
disp(flag);
assert(flag==false, 'CheckMoments did not recognize a valid moment sequence!');
disp('Input:');
disp('------');
M = [1.3, 2.4, 6.03, 20.5, 89.5];
disp('M = ');
disp(M);
disp('Test:');
disp('-----');
disp('flag = CheckMoments(M);:');
flag = CheckMoments(M);
disp('flag = ');
disp(flag);
assert(flag==true, 'CheckMoments did not recognize an invalid moment sequence!');

