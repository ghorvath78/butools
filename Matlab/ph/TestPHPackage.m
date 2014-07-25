format short g;

disp('---BuTools: Continous PH function test file---');

disp('Enable the verbose messages with the BuToolsVerbose flag')
global BuToolsVerbose;
BuToolsVerbose = true

disp('Enable input parameter checking with the BuToolsCheckInput flag')
global BuToolsCheckInput;
BuToolsCheckInput = true


disp('----------------------------------------------------------------------------');
help MomentsFromME

disp('Input:');
disp('------');

a = [0.2, 0.3, 0.5]
A = [-1,0,0;0,-3,2;0,-2,-3]

disp('Test:');
disp('-----');

disp('MomentsFromME(a,A):');
moms=MomentsFromME(a,A);
disp(moms);

assert(CheckMoments(moms)==1, 'The function returned invalid moments!');

disp('Test:');
disp('-----');

disp('MomentsFromME(a,A,9):');
moms=MomentsFromME(a,A,9);
disp(moms);

assert(CheckMoments(moms)==1, 'The function returned invalid moments!');

disp('----------------------------------------------------------------------------');
help MomentsFromPH

disp('Input:');
disp('------');

a=[0.1 0.9 0]
A=[-6.2 2 0; 2 -9 1; 1 0 -3]

disp('Test:');
disp('-----');

disp('MomentsFromPH(a,A,5)');
moms=MomentsFromPH(a,A,5);
disp(moms);

assert(CheckMoments(moms)==1, 'The function returned invalid moments!');

disp('----------------------------------------------------------------------------');
help CdfFromME

disp('Input:');
disp('------');

a = [0.2, 0.3, 0.5]
A = [-1,0,0;0,-3,2;0,-2,-3]

disp('Test:');
disp('-----');

disp('cdf=CdfFromME(a,A,(0:0.01:5)''):');
cdf = CdfFromME(a,A,(0:0.01:5)');

assert(all(diff(cdf)>0), 'The cdf is not increasing monotonously!');
assert(abs(sum(1-cdf)*0.01 - MomentsFromME(a,A,1))<0.01, 'The mean computed from the cdf does not match the theoretical result!');

disp('----------------------------------------------------------------------------');
help CdfFromPH

disp('Input:');
disp('------');

a=[0.1 0.9 0]
A=[-6.2 2 0; 2 -9 1; 1 0 -3]

disp('Test:');
disp('-----');

disp('cdf=CdfFromPH(a,A,(0:0.002:3)''):');
cdf = CdfFromPH(a,A,(0:0.002:3)');

assert(all(diff(cdf)>0), 'The cdf is not increasing monotonously!');
assert(abs(sum(1-cdf)*0.002 - MomentsFromPH(a,A,1))<0.01, 'The mean computed from the cdf does not match the theoretical result!');

disp('----------------------------------------------------------------------------');
help PdfFromME

disp('Input:');
disp('------');

a = [0.2, 0.3, 0.5]
A = [-1,0,0;0,-3,2;0,-2,-3]
x = (0:0.01:5)';

disp('Test:');
disp('-----');

disp('pdf=PdfFromME(a,A,x):');
pdf = PdfFromME(a,A,x);

assert(all(pdf)>=0, 'The pdf is negative!');
assert(abs(pdf'*x*0.01 - MomentsFromME(a,A,1))<0.01, 'The mean computed from the pdf does not match the theoretical result!');

disp('----------------------------------------------------------------------------');
help PdfFromPH

disp('Input:');
disp('------');

a = [0.1 0.9 0]
A = [-6.2 2 0; 2 -9 1; 1 0 -3]
x = (0:0.002:3)';

disp('Test:');
disp('-----');

disp('pdf=PdfFromPH(a,A,x):');
pdf = PdfFromPH(a,A,x);

assert(all(pdf)>=0, 'The pdf is negative!');
assert(abs(pdf'*x*0.002 - MomentsFromPH(a,A,1))<0.002, 'The mean computed from the pdf does not match the theoretical result!');

disp('----------------------------------------------------------------------------');
help IntervalPdfFromPH

disp('Input:');
disp('------');

a = [0.1 0.9 0]
A = [-6.2 2 0; 2 -9 1; 1 0 -3]
x = (0:0.002:3.002)';

disp('Test:');
disp('-----');

disp('[x,y]=IntervalPdfFromPH(a,A,x):');
[x,y] = IntervalPdfFromPH(a,A,x);

assert(all(y)>=0, 'The interval pdf is negative!');
assert(abs(y'*x*0.002 - MomentsFromPH(a,A,1))<0.002, 'The mean computed from the interval pdf does not match the theoretical result!');

disp('----------------------------------------------------------------------------');
help RandomPH

disp('Test:');
disp('-----');
   
disp('[a,A]=RandomPH(3,8,4):');
[a,A]=RandomPH(3,8,4)

assert(CheckPHRepresentation(a,A), 'RandomPH failed to return a valid PH representation!');
assert(max(abs(MomentsFromPH(a,A,1)-8))<1e-14, 'RandomPH failed to match the given mean value!');

disp('----------------------------------------------------------------------------');
help CheckMERepresentation

disp('Input:');
disp('------');

a=[-0.2 0.2]
A=[1 -1; 1 -2]

disp('Test:');
disp('-----');

disp('CheckMERepresentation(a,A):');
flag=CheckMERepresentation(a,A);
disp(flag);

assert(flag==0, 'CheckMERepresentation did not detect that the initial vector is invalid!');

disp('Input:');
disp('------');

a=[-0.2 0.4 0.8]
A=[-2 0 3; 0 -1 1; 0 -1 -1]

disp('Test:');
disp('-----');

disp('CheckMERepresentation(a,A):');
flag=CheckMERepresentation(a,A);
disp(flag);

assert(flag==0, 'CheckMERepresentation did not detect that the dominant eigenvalue is invalid!');

disp('Input:');
disp('------');

a = [0.2, 0.3, 0.5]
A = [-1,0,0;0,-3,2;0,-2,-3]

disp('Test:');
disp('-----');

disp('CheckMERepresentation(a,A):');
flag=CheckMERepresentation(a,A);
disp(flag);

assert(flag==1, 'CheckMERepresentation did not recognize that the given ME representation is valid!');

disp('----------------------------------------------------------------------------');
help CheckPHRepresentation

disp('Input:');
disp('------');

a=[0.2]
A=[-1 1; 1 -2]

disp('Test:');
disp('-----');

disp('CheckPHRepresentation(a,A):');
flag=CheckPHRepresentation(a,A);
disp(flag);

assert(flag==0, 'CheckPHRepresentation did not recognize the wrong input dimensions!');

disp('Input:');
disp('------');

a=[0.2 0.7]
A=[-1 1; 1 -2]

disp('Test:');
disp('-----');

disp('CheckPHRepresentation(a,A):');
flag=CheckPHRepresentation(a,A);
disp(flag);

assert(flag==1, 'CheckPHRepresentation did not recognize that the given PH representation is valid!');

disp('----------------------------------------------------------------------------');
help CheckMEPositiveDensity

disp('Input:');
disp('------');

a = [0.2, 0.3, 0.5]
A = [-1,0,0;0,-3,2;0,-2,-3]

disp('Test:');
disp('-----');

disp('CheckMEPositiveDensity(a,A):');
flag=CheckMEPositiveDensity(a,A);
disp(flag);

assert(flag==1, 'CheckMEPositiveDensity did not recognize that the given ME distribution has positive density!');

disp('Input:');
disp('------');

a = [0.2, 0.3, 0.5]
A = [-1,0,0;0,-3,2.9;0,-2.9,-3]

disp('Test:');
disp('-----');

disp('CheckMEPositiveDensity(a,A):');
flag=CheckMEPositiveDensity(a,A);
disp(flag);

assert(flag==0, 'CheckMEPositiveDensity did not recognize that the given ME distribution does not have positive density!');

disp('----------------------------------------------------------------------------');
help APHFrom3Moments

disp('Input:');
disp('------');

moms = [10,125,8400]
normMoms = NormMomsFromMoms (moms)

disp('Test:');
disp('-----');

disp('[a,A]=APHFrom3Moments(moms(1), moms(2), moms(3))');
disp('size(A,1)');
[a,A]=APHFrom3Moments(moms);
disp(size(A,1));

disp('MomentsFromPH(a,A,3):');
phmoms = MomentsFromPH(a,A,3);
disp(phmoms);

assert(all(abs((phmoms-moms)./moms)<1e-12), 'APHFrom3Moments failed to match the given moments!');

disp('Input:');
disp('------');

moms = [10,525,31400]
normMoms = NormMomsFromMoms (moms)

disp('Test:');
disp('-----');

disp('[a,A]=APHFrom3Moments(moms(1), moms(2), moms(3))');
disp('size(A,1)');
[a,A]=APHFrom3Moments(moms);
disp(size(A,1));

disp('MomentsFromPH(a,A,3):');
phmoms = MomentsFromPH(a,A,3);
disp(phmoms);

assert(all(abs((phmoms-moms)./moms)<1e-12), 'APHFrom3Moments failed to match the given moments!');


disp('----------------------------------------------------------------------------');
help PH2From3Moments

disp('Input:');
disp('------');

moms = [10,160,3500]

disp('Test:');
disp('-----');

disp('[a,A]=PH2From3Moments(moms)');
[a,A]=PH2From3Moments(moms)

disp('MomentsFromPH(a,A,3):');
phmoms = MomentsFromPH(a,A,3);
disp(phmoms);

assert(all(abs((phmoms-moms)./moms)<1e-12), 'PH2From3Moments failed to match the given moments!');

disp('Input:');
disp('------');

moms = [10,260,13500]

disp('Test:');
disp('-----');

disp('[a,A]=PH2From3Moments(moms)');
[a,A]=PH2From3Moments(moms)

disp('MomentsFromPH(a,A,3):');
phmoms = MomentsFromPH(a,A,3);
disp(phmoms);

assert(all(abs((phmoms-moms)./moms)<1e-12), 'PH2From3Moments failed to match the given moments!');

disp('----------------------------------------------------------------------------');
help PH3From5Moments

disp('Input:');
disp('------');

a=[0.1 0.9 0];
A=[-6.2 2 0; 2 -9 1; 1 0 -3];
moms = MomentsFromPH(a,A)

disp('Test:');
disp('-----');

disp('[a,A]=PH3From5Moments(moms)');
[a,A]=PH3From5Moments(moms)

disp('MomentsFromPH(a,A,5):');
phmoms = MomentsFromPH(a,A,5);
disp(phmoms);

assert(all(abs((phmoms-moms)./moms)<1e-12), 'PH3From5Moments failed to match the given moments!');

disp('Input:');
disp('------');

a = [0.2, 0.3, 0.5];
A = [-1,0,0;0,-3,0.5;0,-0.5,-3];
moms = MomentsFromME(a,A)

disp('Test:');
disp('-----');

disp('[a,A]=PH3From5Moments(moms)');
[a,A]=PH3From5Moments(moms)

disp('MomentsFromPH(a,A,5):');
phmoms = MomentsFromPH(a,A,5);
disp(phmoms);

assert(all(abs((phmoms-moms)./moms)<1e-12), 'PH3From5Moments failed to match the given moments!');

disp('----------------------------------------------------------------------------');
help MEFromMoments

disp('Input:');
disp('------');

a=[0.1 0.9 0];
A=[-6.2 2 0; 2 -9 1; 1 0 -3];
moms=MomentsFromPH(a,A,5)

disp('Test:');
disp('-----');

disp('[a,A]=MEFromMoments(moms):');
[a,A]=MEFromMoments(moms)

disp('MomentsFromME(a,A,5):');
memoms = MomentsFromME(a,A,5);
disp(memoms);

assert(all(abs((memoms-moms)./moms)<1e-12), 'MEFromMoments failed to match the given moments!');

disp('----------------------------------------------------------------------------');
help APH2ndMomentLowerBound

disp('Input:');
disp('------');

mean = 1.9
n = 4

disp('Test:');
disp('-----');

disp('mom2 = APH2ndMomentLowerBound(mean,n):');
mom2 = APH2ndMomentLowerBound(mean,n);
disp(mom2);

cv2 = mom2/mean^2-1

assert(abs(cv2-1/n)<1e-14, 'APH2ndMomentLowerBound did not give the expected result!');


disp('----------------------------------------------------------------------------');
help APH3rdMomentLowerBound
help APH3rdMomentUpperBound

disp('Input:');
disp('------');

mean = 1.9
mom2 = 5
n = 3

disp('Test:');
disp('-----');

disp('mom3lower = APH3rdMomentLowerBound(mean,mom2,n):');
mom3lower = APH3rdMomentLowerBound(mean,mom2,n);
disp(mom3lower);
disp('mom3upper = APH3rdMomentUpperBound(mean,mom2,n):');
mom3upper = APH3rdMomentUpperBound(mean,mom2,n);
disp(mom3upper);

assert(mom3upper>mom3lower, 'Lower bound is larger than the upper bound!');

disp('Input:');
disp('------');

mean = 1.9
mom2 = 5
n = 4

disp('Test:');
disp('-----');

disp('mom3lower = APH3rdMomentLowerBound(mean,mom2,n):');
mom3lower = APH3rdMomentLowerBound(mean,mom2,n);
disp(mom3lower);
disp('mom3upper = APH3rdMomentUpperBound(mean,mom2,n):');
mom3upper = APH3rdMomentUpperBound(mean,mom2,n);
disp(mom3upper);

assert(mom3upper>mom3lower, 'Lower bound is larger than the upper bound!');
assert(mom3upper==inf, 'Upper bound must be infinity with 4 phases!');

disp('----------------------------------------------------------------------------');
help CanonicalFromPH2

disp('Input:');
disp('------');

a=[0.12 0.88]
A=[-1.28 0; 3.94 -3.94]

disp('Test:');
disp('-----');

disp('[b,B]=CanonicalFromPH2(a,A):');
[b,B]=CanonicalFromPH2(a,A)

C=SimilarityMatrix(A,B);
err1 = max(max(abs(A*C-C*B)));
err2 = max(abs(a*C-b));

disp('Transformation errors:');
disp(err1);
disp(err2);

assert(err1<1e-12 && err2<1e-12, 'Transformation to canonical PH(2) failed!');

disp('----------------------------------------------------------------------------');
help CanonicalFromPH3

disp('Input:');
disp('------');

a=[0.1 0.9 0]
A=[-6.2 2 0; 2 -9 1; 1 0 -3]

disp('Test:');
disp('-----');

disp('[b,B]=CanonicalFromPH3(a,A):');
[b,B]=CanonicalFromPH3(a,A)

C=SimilarityMatrix(A,B);
err1 = max(max(abs(A*C-C*B)));
err2 = max(abs(a*C-b));

disp('Transformation errors:');
disp(err1);
disp(err2);

assert(err1<1e-12 && err2<1e-12, 'Transformation to canonical PH(3) failed!');

disp('----------------------------------------------------------------------------');
help AcyclicPHFromME

disp('Input:');
disp('------');

a=[-0.4 1.4 0]
A=[-4 1 1; 0 -2 1; 1 0 -8]

disp('Test:');
disp('-----');

disp('[b,B]=AcyclicPHFromME(a,A):');
[b,B]=AcyclicPHFromME(a,A)

disp('Moments of (a,A) and (b,B)');
ma=MomentsFromME(a,A,5);
mb=MomentsFromME(b,B,5);
disp(ma);
disp(mb);

assert(norm((ma-mb)./ma)<1e-7, 'Transformation to acyclic representation failed!');

disp('----------------------------------------------------------------------------');
help MonocyclicPHFromME

disp('Input:');
disp('------');

a = [0.2, 0.3, 0.5]
A = [-1,0,0;0,-3,2;0,-2,-3]

disp('Test:');
disp('-----');

disp('[b,B]=MonocyclicPHFromME(a,A):');
[b,B]=MonocyclicPHFromME(a,A)

disp('Moments of (a,A) and (b,B)');
ma=MomentsFromME(a,A,5);
mb=MomentsFromME(b,B,5);
disp(ma);
disp(mb);

assert(norm((ma-mb)./ma)<1e-7, 'Transformation to monocyclic representation failed!');

disp('----------------------------------------------------------------------------');
help PHFromME

disp('Input:');
disp('------');

a=[-0.4 1.4]
A=[-3.8 2; 2 -9]

disp('Test:');
disp('-----');

disp('CheckMERepresentation(a,A):');
flag=CheckMERepresentation(a,A);
disp(flag);

disp('CheckPHRepresentation(a,A):');
flag=CheckPHRepresentation(a,A);
disp(flag);

disp('[b,B]=PHFromME(a,A)');
[b,B]=PHFromME(a,A)

disp('Check the obtained PH, CheckPHRepresentation(b,B):');
flag=CheckPHRepresentation(b,B);
disp(flag);

C=SimilarityMatrix(A,B);
err1 = max(max(abs(A*C-C*B)));
err2 = max(abs(a*C-b));

assert(flag && err1<1e-12 && err2<1e-12, 'Transformation to PH failed!');

disp('Input:');
disp('------');

a=[-0.5 1.5]
A=[-3.8 2; 2 -9]

disp('Test:');
disp('-----');

disp('CheckMERepresentation(a,A):');
flag=CheckMERepresentation(a,A);
disp(flag);

disp('CheckPHRepresentation(a,A):');
flag=CheckPHRepresentation(a,A);
disp(flag);

disp('[b,B]=PHFromME(a,A)');
[b,B]=PHFromME(a,A)

disp('Check the obtained PH, CheckPHRepresentation(b,B):');
flag=CheckPHRepresentation(b,B);
disp(flag);

C=SimilarityMatrix(A,B);
err1 = max(max(abs(A*C-C*B)));
err2 = max(abs(a*C-b));

assert(flag && err1<1e-12 && err2<1e-12, 'Transformation to PH failed!');

disp('----------------------------------------------------------------------------');
help MEOrder

disp('Input:');
disp('------');

a=[1 1 1 1 1 1]/6
A=[-1 0 0 0 0 0; 0.5 -2 1 0 0 0; 1 0 -3 1 0 0; 1 0 1 -4 1 0; ...
    4 0 0 0 -5 0; 5 0 0 0 0 -6]

disp('Test:');
disp('-----');

disp('co=MEOrder(a,A,''cont'')');
co=MEOrder(a,A,'cont');
disp(co);
disp('oo=MEOrder(a,A,''obs'')');
oo=MEOrder(a,A,'obs');
disp(oo);
disp('coo=MEOrder(a,A,''obscont'')');
coo=MEOrder(a,A,'obscont');
disp(coo);
disp('coo=MEOrder(a,A,''moment'')');
mo=MEOrder(a,A,'moment');
disp(mo);

assert(co==2, 'Wrong controllability order returned!');
assert(oo==6, 'Wrong observability order returned!');
assert(coo==2, 'The minimum of the controllability and observability order is wrong!');
assert(mo==2, 'Wrong moment order returned!');

disp('Input:');
disp('------');

a=[2 1]/3
A=[-1 1; 0 -3]

disp('Test:');
disp('-----');

disp('co=MEOrder(a,A,''cont'')');
co=MEOrder(a,A,'cont');
disp(co);
disp('oo=MEOrder(a,A,''obs'')');
oo=MEOrder(a,A,'obs');
disp(oo);
disp('coo=MEOrder(a,A,''obscont'')');
coo=MEOrder(a,A,'obscont');
disp(coo);
disp('coo=MEOrder(a,A,''moment'')');
mo=MEOrder(a,A,'moment');
disp(mo);

assert(co==2, 'Wrong controllability order returned!');
assert(oo==1, 'Wrong observability order returned!');
assert(coo==1, 'The minimum of the controllability and observability order is wrong!');
assert(mo==1, 'Wrong moment order returned!');

disp('Input:');
disp('------');

b = [0.2, 0.3, 0.5];
B = [-1,0,0;0,-3,1;0,-1,-3];
[a,A] = MonocyclicPHFromME(b,B)

disp('Test:');
disp('-----');

disp('co=MEOrder(a,A,''cont'')');
co=MEOrder(a,A,'cont');
disp(co);
disp('oo=MEOrder(a,A,''obs'')');
oo=MEOrder(a,A,'obs');
disp(oo);
disp('coo=MEOrder(a,A,''obscont'')');
coo=MEOrder(a,A,'obscont');
disp(coo);
disp('coo=MEOrder(a,A,''moment'')');
mo=MEOrder(a,A,'moment');
disp(mo);

assert(co==9, 'Wrong controllability order returned!');
assert(oo==3, 'Wrong observability order returned!');
assert(coo==3, 'The minimum of the controllability and observability order is wrong!');
assert(mo==3, 'Wrong moment order returned!');

disp('----------------------------------------------------------------------------');
help MEOrderFromMoments

disp('Input:');
disp('------');

a=[0.1 0.9 0]
A=[-6.2 2 0; 2 -9 1; 1 0 -3]

disp('Test:');
disp('-----');

disp('moms=MomentsFromME(a,A)');
disp('mo=MEOrderFromMoments(moms):');
moms=MomentsFromME(a,A);
mo = MEOrderFromMoments(moms);
disp(mo);

assert(mo==3, 'Wrong moment order returned!');

disp('Input:');
disp('------');

b = [0.2, 0.3, 0.5];
B = [-1,0,0;0,-3,2;0,-2,-3];
[a,A] = MonocyclicPHFromME(b,B)

disp('Test:');
disp('-----');

disp('moms=MomentsFromME(a,A)');
disp('mo=MEOrderFromMoments(moms):');
moms=MomentsFromME(a,A);
mo = MEOrderFromMoments(moms);
disp(mo);

assert(mo==3, 'Wrong moment order returned!');

disp('----------------------------------------------------------------------------');
help MinimalRepFromME

disp('Input:');
disp('------');

a=[1 1 1 1 1 1]/6
A=[-1 0 0 0 0 0; 0.5 -2 1 0 0 0; 1 0 -3 1 0 0; 1 0 1 -4 1 0; ...
    4 0 0 0 -5 0; 5 0 0 0 0 -6]

disp('Test:');
disp('-----');

disp('[b,B]=MinimalRepFromME(a,A,''cont'')');
[b,B]=MinimalRepFromME(a,A,'cont');
disp(b);
disp(B);

assert(length(b)==2, 'Non-minimal representation returned based on controllability!');

disp('[b,B]=MinimalRepFromME(a,A,''obs'')');
[b,B]=MinimalRepFromME(a,A,'obs');
disp(b);
disp(B);

assert(length(b)==6, 'Non-minimal representation returned based on observability!');

disp('[b,B]=MinimalRepFromME(a,A,''obscont'')');
[b,B]=MinimalRepFromME(a,A,'obscont');
disp(b);
disp(B);

assert(length(b)==2, 'Non-minimal representation returned based on observability and controllability!');

disp('[b,B]=MinimalRepFromME(a,A,''moment'')');
[b,B]=MinimalRepFromME(a,A,'moment');
disp(b);
disp(B);

assert(length(b)==2, 'Non-minimal representation returned based on the moments!');

disp('Input:');
disp('------');

a=[2 1]/3
A=[-1 1; 0 -3]

disp('Test:');
disp('-----');

disp('[b,B]=MinimalRepFromME(a,A,''cont'')');
[b,B]=MinimalRepFromME(a,A,'cont');
disp(b);
disp(B);

assert(length(b)==2, 'Non-minimal representation returned based on controllability!');

disp('[b,B]=MinimalRepFromME(a,A,''obs'')');
[b,B]=MinimalRepFromME(a,A,'obs');
disp(b);
disp(B);

assert(length(b)==1, 'Non-minimal representation returned based on observability!');

disp('[b,B]=MinimalRepFromME(a,A,''obscont'')');
[b,B]=MinimalRepFromME(a,A,'obscont');
disp(b);
disp(B);

assert(length(b)==1, 'Non-minimal representation returned based on observability and controllability!');

disp('[b,B]=MinimalRepFromME(a,A,''moment'')');
[b,B]=MinimalRepFromME(a,A,'moment');
disp(b);
disp(B);

assert(length(b)==1, 'Non-minimal representation returned based on the moments!');

disp('Input:');
disp('------');

b = [0.2, 0.3, 0.5];
B = [-1,0,0;0,-3,1;0,-1,-3];
[a,A] = MonocyclicPHFromME(b,B)

disp('[b,B]=MinimalRepFromME(a,A,''cont'')');
[b,B]=MinimalRepFromME(a,A,'cont');
disp(b);
disp(B);

assert(length(b)==length(a), 'Non-minimal representation returned based on controllability!');

disp('[b,B]=MinimalRepFromME(a,A,''obs'')');
[b,B]=MinimalRepFromME(a,A,'obs');
disp(b);
disp(B);
% check similarity:
C = SimilarityMatrix(B,A);
sim = norm(B*C-C*A) + norm(b*C-a);

assert(length(b)==3 && sim<1e-14, 'Non-minimal representation returned based on observability!');

disp('[b,B]=MinimalRepFromME(a,A,''obscont'')');
[b,B]=MinimalRepFromME(a,A,'obscont');
disp(b);
disp(B);
% check similarity:
C = SimilarityMatrix(B,A);
sim = norm(B*C-C*A) + norm(b*C-a);

assert(length(b)==3 && sim<1e-14, 'Non-minimal representation returned based on observability and controllability!');

disp('[b,B]=MinimalRepFromME(a,A,''moment'')');
[b,B]=MinimalRepFromME(a,A,'moment');
disp(b);
disp(B);
% check similarity:
C = SimilarityMatrix(B,A);
sim = norm(B*C-C*A) + norm(b*C-a);

assert(length(b)==3 && sim<1e-14, 'Non-minimal representation returned based on the moments!');

disp('----------------------------------------------------------------------------');
help SamplesFromPH

disp('Input:');
disp('------');

a=[0.1 0.9 0]
A=[-6.2 2 0; 2 -9 1; 1 0 -3]

disp('Test:');
disp('-----');

disp('x=SamplesFromPH(a,A,1000)');
x=SamplesFromPH(a,A,1000);

disp('Moments from the samples:');
mt = MarginalMomentsFromTrace(x,3);
disp(mt);

disp('Moments from the PH:');
mp = MomentsFromPH(a,A,3);
disp(mp);
