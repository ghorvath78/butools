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
format compact
delete('/home/gabor/github/butools/test/docex/Trace_matlab.docex');
diary('/home/gabor/github/butools/test/docex/Trace_matlab.docex');
disp('=== CdfFromTrace ===')
disp('>>> D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];');
D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];
disp('>>> D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];');
D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];
disp('>>> tr = SamplesFromMAP(D0, D1, 1000000);');
tr = SamplesFromMAP(D0, D1, 1000000);
disp('>>> [x, y] = CdfFromTrace(tr);');
[x, y] = CdfFromTrace(tr);
disp('>>> plot(x, y)');
meandiff = abs(dot(diff(x), 1.-y(1:end-1))-mean(tr))/mean(tr);
disp('=== PdfFromTrace ===')
disp('>>> D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];');
D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];
disp('>>> D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];');
D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];
disp('>>> x = (0.0:0.01:0.5);');
x = (0.0:0.01:0.5);
disp('>>> tr = SamplesFromMAP(D0, D1, 1000000);');
tr = SamplesFromMAP(D0, D1, 1000000);
disp('>>> [x, y] = PdfFromTrace(tr, x);');
[x, y] = PdfFromTrace(tr, x);
disp('>>> [a, A] = MarginalDistributionFromMAP(D0, D1);');
[a, A] = MarginalDistributionFromMAP(D0, D1);
disp('>>> [xm, ym] = IntervalPdfFromPH(a, A, x);');
[xm, ym] = IntervalPdfFromPH(a, A, x);
disp('>>> plot(x, y, xm, ym)');
meandiff = abs(dot(x(1:end-1), (diff(x).*y(1:end-1)))-mean(tr))/mean(tr);
disp('=== MarginalMomentsFromTrace ===')
disp('>>> D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];');
D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];
disp('>>> D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];');
D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];
disp('>>> tr = SamplesFromMAP(D0, D1, 1000000);');
tr = SamplesFromMAP(D0, D1, 1000000);
disp('>>> moms = MarginalMomentsFromTrace(tr, 3);');
moms = MarginalMomentsFromTrace(tr, 3);
disp('>>> disp(moms);');
disp(moms);
disp('>>> mmoms = MarginalMomentsFromMAP(D0, D1, 3);');
mmoms = MarginalMomentsFromMAP(D0, D1, 3);
disp('>>> disp(mmoms);');
disp(mmoms);
disp('=== LagCorrelationsFromTrace ===')
disp('>>> D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];');
D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];
disp('>>> D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];');
D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];
disp('>>> tr = SamplesFromMAP(D0, D1, 1000000);');
tr = SamplesFromMAP(D0, D1, 1000000);
disp('>>> acf = LagCorrelationsFromTrace(tr, 10);');
acf = LagCorrelationsFromTrace(tr, 10);
disp('>>> disp(acf);');
disp(acf);
disp('>>> macf = LagCorrelationsFromMAP(D0, D1, 10);');
macf = LagCorrelationsFromMAP(D0, D1, 10);
disp('>>> disp(macf);');
disp(macf);
disp('=== LagkJointMomentsFromTrace ===')
disp('>>> D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];');
D0 = [-18., 1., 4.; 2., -18., 7.; 1., 3., -32.];
disp('>>> D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];');
D1 = [12., 1., 0.; 1., 8., 0.; 2., 1., 25.];
disp('>>> tr = SamplesFromMAP(D0, D1, 1000000);');
tr = SamplesFromMAP(D0, D1, 1000000);
disp('>>> Nm1 = LagkJointMomentsFromTrace(tr, 3, 1);');
Nm1 = LagkJointMomentsFromTrace(tr, 3, 1);
disp('>>> disp(Nm1);');
disp(Nm1);
disp('>>> mNm1 = LagkJointMomentsFromMAP(D0, D1, 3, 1);');
mNm1 = LagkJointMomentsFromMAP(D0, D1, 3, 1);
disp('>>> disp(mNm1);');
disp(mNm1);
disp('>>> Nm2 = LagkJointMomentsFromTrace(tr, 3, 2);');
Nm2 = LagkJointMomentsFromTrace(tr, 3, 2);
disp('>>> disp(Nm2);');
disp(Nm2);
disp('>>> mNm2 = LagkJointMomentsFromMAP(D0, D1, 3, 2);');
mNm2 = LagkJointMomentsFromMAP(D0, D1, 3, 2);
disp('>>> disp(mNm2);');
disp(mNm2);
disp('=== CdfFromWeightedTrace ===')
disp('>>> wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34];');
wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34];
disp('>>> weights = [12., 1., 34., 23., 8., 2.];');
weights = [12., 1., 34., 23., 8., 2.];
disp('>>> [x, y] = CdfFromWeightedTrace(wtrace, weights);');
[x, y] = CdfFromWeightedTrace(wtrace, weights);
disp('>>> plot(x, y)');
disp('=== PdfFromWeightedTrace ===')
disp('>>> wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34];');
wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34];
disp('>>> weights = [12., 1., 34., 23., 8., 2.];');
weights = [12., 1., 34., 23., 8., 2.];
disp('>>> x = (0.0:0.1:3.0);');
x = (0.0:0.1:3.0);
disp('>>> [x, y] = PdfFromWeightedTrace(wtrace, weights, x);');
[x, y] = PdfFromWeightedTrace(wtrace, weights, x);
disp('>>> plot(x, y)');
disp('=== MarginalMomentsFromWeightedTrace ===')
disp('>>> wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34];');
wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34];
disp('>>> weights = [12., 1., 34., 23., 8., 2.];');
weights = [12., 1., 34., 23., 8., 2.];
disp('>>> moms = MarginalMomentsFromWeightedTrace(wtrace, weights, 3);');
moms = MarginalMomentsFromWeightedTrace(wtrace, weights, 3);
disp('>>> disp(moms);');
disp(moms);

