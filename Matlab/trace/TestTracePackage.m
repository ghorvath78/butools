format short g;

disp('---BuTools: Trace processing function test file---');

disp('Enable the verbose messages with the BuToolsVerbose flag')
global BuToolsVerbose;
BuToolsVerbose = true

disp('Enable input parameter checking with the BuToolsCheckInput flag')
global BuToolsCheckInput;
BuToolsCheckInput = true

disp('----------------------------------------------------------------------------');
disp('Generating trace file for the tests...');

D0 = [-18 1 4; 2 -18 7; 1 3 -32];
D1 = [12 1 0; 1 8 0; 2 1 25]; 
tr = SamplesFromMAP(D0,D1,1000000);
mmoms = MarginalMomentsFromMAP(D0,D1,5);
macf = LagCorrelationsFromMAP(D0,D1,10);
mNm1 = LagkJointMomentsFromMAP(D0,D1,3,1);
mNm2 = LagkJointMomentsFromMAP(D0,D1,3,2);
[a,A]=MarginalDistributionFromMAP(D0,D1);

disp('----------------------------------------------------------------------------');
help CdfFromTrace

disp('Test:');
disp('-----');

disp('[x,y]=CdfFromTrace(tr):');
[x,y]=CdfFromTrace(tr);
plot(x,y);

meandiff = abs(diff(x)'*(1.-y(1:end-1)) - mean(tr))/mean(tr);
assert(all(diff(y)>=0) && y(end)<=1 && y(1)>=0, 'CdfFromTrace returned a wrong cdf!');
assert(meandiff<1e-2, 'The mean obtained from the cdf returned by CdfFromTrace does not match the trace mean!');

disp('----------------------------------------------------------------------------');
help PdfFromTrace

disp('Test:');
disp('-----');

disp('[x,y]=PdfFromTrace(tr, 0:0.01:0.5):');
[x,y]=PdfFromTrace(tr, 0:0.01:0.5);
[xm,ym]=IntervalPdfFromPH(a, A, 0:0.01:0.5);
plot(x,y,xm,ym);

meandiff = abs(x(1:end-1)'*(diff(x).*y(1:end-1)) - mean(tr))/mean(tr);

assert(all(y>=0), 'PdfFromTrace returned a wrong pdf!');
assert(meandiff<1e-2, 'The mean obtained from the pdf returned by PdfFromTrace does not match the trace mean!');

disp('----------------------------------------------------------------------------');
help MarginalMomentsFromTrace

disp('Test:');
disp('-----');

disp('moms=MarginalMomentsFromTrace(tr, 3):');
moms=MarginalMomentsFromTrace(tr, 3);
moms
mmoms(1:3)

momdiff = norm((mmoms(1:3)-moms)./mmoms(1:3));

assert(momdiff<1e-1, 'Moments from MarginalMomentsFromTrace are far from the theoretical moments of the trace!');

disp('----------------------------------------------------------------------------');
help LagCorrelationsFromTrace

disp('Test:');
disp('-----');

disp('acf=LagCorrelationsFromTrace(tr, 10):');
acf=LagCorrelationsFromTrace(tr, 10);
acf
macf
plot([acf,macf]);

acfdiff = norm(acf-macf);

assert(acfdiff<1e-1, 'Autocorrelations from LagCorrelationsFromTrace are far from the theoretical autocorrelations of the trace!');

disp('----------------------------------------------------------------------------');
help LagkJointMomentsFromTrace


disp('Test:');
disp('-----');

disp('Nm=LagkJointMomentsFromTrace(tr,3,1):');
Nm=LagkJointMomentsFromTrace(tr,3,1)
mNm1

Nmdiff = norm((Nm-mNm1)./mNm1);
assert(Nmdiff<1e-1, 'Moments from LagkJointMomentsFromTrace are far from the theoretical moments of the trace!');

disp('Test:');
disp('-----');

disp('Nm=LagkJointMomentsFromTrace(tr,3,2):');
Nm=LagkJointMomentsFromTrace(tr,3,2)
mNm2

Nmdiff = norm((Nm-mNm2)./mNm2);
assert(Nmdiff<1e-1, 'Moments from LagkJointMomentsFromTrace are far from the theoretical moments of the trace!');

disp('----------------------------------------------------------------------------');
help CdfFromWeightedTrace

disp('Input:');
disp('------');
wtr=[0.12; 1.23; 0.546; 0.6765; 1.34; 2.34];
wei=[12; 1; 34; 23; 8; 2];

disp('Test:');
disp('-----');

disp('[x,y]=CdfFromWeightedTrace(wtr,wei):');
[x,y]=CdfFromWeightedTrace(wtr,wei);
plot(x,y);

assert(all(diff(y)>=0) && y(end)<=1 && y(1)>=0, 'CdfFromWeightedTrace returned a wrong cdf!');

disp('----------------------------------------------------------------------------');
help PdfFromWeightedTrace

disp('Test:');
disp('-----');

disp('[x,y]=PdfFromWeightedTrace(wtr, wei 0:0.1:3):');
[x,y]=PdfFromWeightedTrace(wtr, wei, 0:0.1:3);
plot(x,y);

assert(all(y>=0), 'PdfFromWeightedTrace returned a wrong pdf!');

disp('----------------------------------------------------------------------------');
help MarginalMomentsFromWeightedTrace

disp('Test:');
disp('-----');

disp('moms=MarginalMomentsFromWeightedTrace(wtr, wei, 3):');
moms=MarginalMomentsFromWeightedTrace(wtr, wei, 3);
moms

