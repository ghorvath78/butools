clear all
run('/home/gabor/github/butools/Matlab/BuToolsInit.m')

disp('---BuTools: Trace package test file---');
disp('Enable the verbose messages with the BuToolsVerbose flag');
global BuToolsVerbose;
BuToolsVerbose = true;
disp('Enable input parameter checking with the BuToolsCheckInput flag');
global BuToolsCheckInput;
BuToolsCheckInput = true;
global BuToolsCheckPrecision;
format short g;
disp('========================================')
disp('Testing BuTools function CdfFromTrace')
disp('Input:');
disp('------');
D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];
disp('D0 = ');
disp(D0);
D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];
disp('D1 = ');
disp(D1);
disp('tr = SamplesFromMAP(D0, D1, 1000000);:');
tr = SamplesFromMAP(D0, D1, 1000000);
disp('Test:');
disp('-----');
disp('[x, y] = CdfFromTrace(tr);:');
[x, y] = CdfFromTrace(tr);
plot(x, y)
meandiff = abs(dot(diff(x), 1.-y(1:end-1))-mean(tr))/mean(tr);
assert(all(diff(y)>=0)&&y(end)<=1&&y(1)>=0, 'CdfFromTrace returned a wrong cdf!');
assert(meandiff<10^-2, 'The mean obtained from the cdf returned by CdfFromTrace does not match the trace mean!');
disp('========================================')
disp('Testing BuTools function PdfFromTrace')
disp('Input:');
disp('------');
D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];
disp('D0 = ');
disp(D0);
D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];
disp('D1 = ');
disp(D1);
x = (0.0:0.01:0.5);
disp('tr = SamplesFromMAP(D0, D1, 1000000);:');
tr = SamplesFromMAP(D0, D1, 1000000);
disp('Test:');
disp('-----');
disp('[x, y] = PdfFromTrace(tr, x);:');
[x, y] = PdfFromTrace(tr, x);
disp('[a, A] = MarginalDistributionFromMAP(D0, D1);:');
[a, A] = MarginalDistributionFromMAP(D0, D1);
disp('[xm, ym] = IntervalPdfFromPH(a, A, x);:');
[xm, ym] = IntervalPdfFromPH(a, A, x);
plot(x, y, xm, ym)
meandiff = abs(dot(x(1:end-1), (diff(x).*y(1:end-1)))-mean(tr))/mean(tr);
assert(all(y>=0), 'PdfFromTrace returned a wrong pdf!');
assert(meandiff<10^-2, 'The mean obtained from the pdf returned by PdfFromTrace does not match the trace mean!');
disp('========================================')
disp('Testing BuTools function MarginalMomentsFromTrace')
disp('Input:');
disp('------');
D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];
disp('D0 = ');
disp(D0);
D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];
disp('D1 = ');
disp(D1);
disp('tr = SamplesFromMAP(D0, D1, 1000000);:');
tr = SamplesFromMAP(D0, D1, 1000000);
disp('Test:');
disp('-----');
disp('moms = MarginalMomentsFromTrace(tr, 3);:');
moms = MarginalMomentsFromTrace(tr, 3);
disp('moms = ');
disp(moms);
disp('mmoms = MarginalMomentsFromMAP(D0, D1, 3);:');
mmoms = MarginalMomentsFromMAP(D0, D1, 3);
disp('mmoms = ');
disp(mmoms);
assert(norm((moms-mmoms)./mmoms)<10^-1, 'Moments from MarginalMomentsFromTrace are far from the theoretical moments of the trace!');
disp('========================================')
disp('Testing BuTools function LagCorrelationsFromTrace')
disp('Input:');
disp('------');
D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];
disp('D0 = ');
disp(D0);
D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];
disp('D1 = ');
disp(D1);
disp('tr = SamplesFromMAP(D0, D1, 1000000);:');
tr = SamplesFromMAP(D0, D1, 1000000);
disp('Test:');
disp('-----');
disp('acf = LagCorrelationsFromTrace(tr, 10);:');
acf = LagCorrelationsFromTrace(tr, 10);
disp('acf = ');
disp(acf);
disp('macf = LagCorrelationsFromMAP(D0, D1, 10);:');
macf = LagCorrelationsFromMAP(D0, D1, 10);
disp('macf = ');
disp(macf);
assert(norm(acf-macf)<10^-1, 'Autocorrelations from LagCorrelationsFromTrace are far from the theoretical autocorrelations of the trace!');
disp('========================================')
disp('Testing BuTools function LagkJointMomentsFromTrace')
disp('Input:');
disp('------');
D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];
disp('D0 = ');
disp(D0);
D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];
disp('D1 = ');
disp(D1);
disp('tr = SamplesFromMAP(D0, D1, 1000000);:');
tr = SamplesFromMAP(D0, D1, 1000000);
disp('Test:');
disp('-----');
disp('Nm1 = LagkJointMomentsFromTrace(tr, 3, 1);:');
Nm1 = LagkJointMomentsFromTrace(tr, 3, 1);
disp('Nm1 = ');
disp(Nm1);
disp('mNm1 = LagkJointMomentsFromMAP(D0, D1, 3, 1);:');
mNm1 = LagkJointMomentsFromMAP(D0, D1, 3, 1);
disp('mNm1 = ');
disp(mNm1);
assert(norm((Nm1-mNm1)./mNm1)<10^-1, 'Moments from LagkJointMomentsFromTrace are far from the theoretical moments of the trace!');
disp('Test:');
disp('-----');
disp('Nm2 = LagkJointMomentsFromTrace(tr, 3, 2);:');
Nm2 = LagkJointMomentsFromTrace(tr, 3, 2);
disp('Nm2 = ');
disp(Nm2);
disp('mNm2 = LagkJointMomentsFromMAP(D0, D1, 3, 2);:');
mNm2 = LagkJointMomentsFromMAP(D0, D1, 3, 2);
disp('mNm2 = ');
disp(mNm2);
assert(norm((Nm2-mNm2)./mNm2)<10^-1, 'Moments from LagkJointMomentsFromTrace are far from the theoretical moments of the trace!');
disp('========================================')
disp('Testing BuTools function CdfFromWeightedTrace')
disp('Input:');
disp('------');
wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34];
disp('wtrace = ');
disp(wtrace);
weights = [12., 1., 34., 23., 8., 2.];
disp('weights = ');
disp(weights);
disp('Test:');
disp('-----');
disp('[x, y] = CdfFromWeightedTrace(wtrace, weights);:');
[x, y] = CdfFromWeightedTrace(wtrace, weights);
plot(x, y)
assert(all(diff(y)>=0)&&y(end)<=1&&y(1)>=0, 'CdfFromWeightedTrace returned a wrong cdf!');
disp('========================================')
disp('Testing BuTools function PdfFromWeightedTrace')
disp('Input:');
disp('------');
wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34];
disp('wtrace = ');
disp(wtrace);
weights = [12., 1., 34., 23., 8., 2.];
disp('weights = ');
disp(weights);
x = (0.0:0.1:3.0);
disp('Test:');
disp('-----');
disp('[x, y] = PdfFromWeightedTrace(wtrace, weights, x);:');
[x, y] = PdfFromWeightedTrace(wtrace, weights, x);
plot(x, y)
assert(all(y>=0), 'PdfFromWeightedTrace returned a wrong pdf!');
disp('========================================')
disp('Testing BuTools function MarginalMomentsFromWeightedTrace')
disp('Input:');
disp('------');
wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34];
disp('wtrace = ');
disp(wtrace);
weights = [12., 1., 34., 23., 8., 2.];
disp('weights = ');
disp(weights);
disp('Test:');
disp('-----');
disp('moms = MarginalMomentsFromWeightedTrace(wtrace, weights, 3);:');
moms = MarginalMomentsFromWeightedTrace(wtrace, weights, 3);
disp('moms = ');
disp(moms);

