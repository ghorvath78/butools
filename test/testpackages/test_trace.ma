ClearAll["Global`*"]
AppendTo[$Path,"/home/gabor/github/butools/Mathematica"];
<<BuTools`
Print["---BuTools: Trace package test file---"//OutputForm];
Print["Enable the verbose messages with the BuToolsVerbose flag"//OutputForm];
BuTools`Verbose = True;
Print["Enable input parameter checking with the BuToolsCheckInput flag"//OutputForm];
BuTools`CheckInput = true;
On[Assert];
Print["========================================"]
Print["Testing BuTools function CdfFromTrace"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};
Print["D0 = "//OutputForm];
Print[D0];
D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};
Print["D1 = "//OutputForm];
Print[D1];
Print["tr = SamplesFromMAP[D0, D1, 1000000];:"//OutputForm];
tr = SamplesFromMAP[D0, D1, 1000000];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["{x, y} = CdfFromTrace[tr];:"//OutputForm];
{x, y} = CdfFromTrace[tr];
ListLinePlot[{Transpose[{x, y}]}]
meandiff = Abs[Dot[Differences[x], 1.-y[[1;;-2]]]-Mean[tr]]/Mean[tr];
Assert[AllTrue[Differences[y],#>=-10^-14 &]&&y[[-1]]<=1&&y[[1]]>=0, "CdfFromTrace returned a wrong cdf!"];
Assert[meandiff<10^-2, "The mean obtained from the cdf returned by CdfFromTrace does not match the trace mean!"];
Print["========================================"]
Print["Testing BuTools function PdfFromTrace"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};
Print["D0 = "//OutputForm];
Print[D0];
D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};
Print["D1 = "//OutputForm];
Print[D1];
x = Range[0.0,0.5,0.01];
Print["tr = SamplesFromMAP[D0, D1, 1000000];:"//OutputForm];
tr = SamplesFromMAP[D0, D1, 1000000];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["{x, y} = PdfFromTrace[tr, x];:"//OutputForm];
{x, y} = PdfFromTrace[tr, x];
Print["{a, A} = MarginalDistributionFromMAP[D0, D1];:"//OutputForm];
{a, A} = MarginalDistributionFromMAP[D0, D1];
Print["{xm, ym} = IntervalPdfFromPH[a, A, x];:"//OutputForm];
{xm, ym} = IntervalPdfFromPH[a, A, x];
ListLinePlot[{Transpose[{x, y}],Transpose[{xm, ym}]}]
meandiff = Abs[Dot[x[[1;;-2]], (Differences[x]*y[[1;;-2]])]-Mean[tr]]/Mean[tr];
Assert[AllTrue[y,#>=-10^-14 &], "PdfFromTrace returned a wrong pdf!"];
Assert[meandiff<10^-2, "The mean obtained from the pdf returned by PdfFromTrace does not match the trace mean!"];
Print["========================================"]
Print["Testing BuTools function MarginalMomentsFromTrace"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};
Print["D0 = "//OutputForm];
Print[D0];
D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};
Print["D1 = "//OutputForm];
Print[D1];
Print["tr = SamplesFromMAP[D0, D1, 1000000];:"//OutputForm];
tr = SamplesFromMAP[D0, D1, 1000000];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["moms = MarginalMomentsFromTrace[tr, 3];:"//OutputForm];
moms = MarginalMomentsFromTrace[tr, 3];
Print["moms = "//OutputForm];
Print[moms];
Print["mmoms = MarginalMomentsFromMAP[D0, D1, 3];:"//OutputForm];
mmoms = MarginalMomentsFromMAP[D0, D1, 3];
Print["mmoms = "//OutputForm];
Print[mmoms];
Assert[Norm[(moms-mmoms)/mmoms]<10^-1, "Moments from MarginalMomentsFromTrace are far from the theoretical moments of the trace!"];
Print["========================================"]
Print["Testing BuTools function LagCorrelationsFromTrace"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};
Print["D0 = "//OutputForm];
Print[D0];
D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};
Print["D1 = "//OutputForm];
Print[D1];
Print["tr = SamplesFromMAP[D0, D1, 1000000];:"//OutputForm];
tr = SamplesFromMAP[D0, D1, 1000000];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["acf = LagCorrelationsFromTrace[tr, 10];:"//OutputForm];
acf = LagCorrelationsFromTrace[tr, 10];
Print["acf = "//OutputForm];
Print[acf];
Print["macf = LagCorrelationsFromMAP[D0, D1, 10];:"//OutputForm];
macf = LagCorrelationsFromMAP[D0, D1, 10];
Print["macf = "//OutputForm];
Print[macf];
Assert[Norm[acf-macf]<10^-1, "Autocorrelations from LagCorrelationsFromTrace are far from the theoretical autocorrelations of the trace!"];
Print["========================================"]
Print["Testing BuTools function LagkJointMomentsFromTrace"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};
Print["D0 = "//OutputForm];
Print[D0];
D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};
Print["D1 = "//OutputForm];
Print[D1];
Print["tr = SamplesFromMAP[D0, D1, 1000000];:"//OutputForm];
tr = SamplesFromMAP[D0, D1, 1000000];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["Nm1 = LagkJointMomentsFromTrace[tr, 3, 1];:"//OutputForm];
Nm1 = LagkJointMomentsFromTrace[tr, 3, 1];
Print["Nm1 = "//OutputForm];
Print[Nm1];
Print["mNm1 = LagkJointMomentsFromMAP[D0, D1, 3, 1];:"//OutputForm];
mNm1 = LagkJointMomentsFromMAP[D0, D1, 3, 1];
Print["mNm1 = "//OutputForm];
Print[mNm1];
Assert[Norm[(Nm1-mNm1)/mNm1]<10^-1, "Moments from LagkJointMomentsFromTrace are far from the theoretical moments of the trace!"];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["Nm2 = LagkJointMomentsFromTrace[tr, 3, 2];:"//OutputForm];
Nm2 = LagkJointMomentsFromTrace[tr, 3, 2];
Print["Nm2 = "//OutputForm];
Print[Nm2];
Print["mNm2 = LagkJointMomentsFromMAP[D0, D1, 3, 2];:"//OutputForm];
mNm2 = LagkJointMomentsFromMAP[D0, D1, 3, 2];
Print["mNm2 = "//OutputForm];
Print[mNm2];
Assert[Norm[(Nm2-mNm2)/mNm2]<10^-1, "Moments from LagkJointMomentsFromTrace are far from the theoretical moments of the trace!"];
Print["========================================"]
Print["Testing BuTools function CdfFromWeightedTrace"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
wtrace = {0.12, 1.23, 0.546, 0.6765, 1.34, 2.34};
Print["wtrace = "//OutputForm];
Print[wtrace];
weights = {12., 1., 34., 23., 8., 2.};
Print["weights = "//OutputForm];
Print[weights];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["{x, y} = CdfFromWeightedTrace[wtrace, weights];:"//OutputForm];
{x, y} = CdfFromWeightedTrace[wtrace, weights];
ListLinePlot[{Transpose[{x, y}]}]
Assert[AllTrue[Differences[y],#>=-10^-14 &]&&y[[-1]]<=1&&y[[1]]>=0, "CdfFromWeightedTrace returned a wrong cdf!"];
Print["========================================"]
Print["Testing BuTools function PdfFromWeightedTrace"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
wtrace = {0.12, 1.23, 0.546, 0.6765, 1.34, 2.34};
Print["wtrace = "//OutputForm];
Print[wtrace];
weights = {12., 1., 34., 23., 8., 2.};
Print["weights = "//OutputForm];
Print[weights];
x = Range[0.0,3.0,0.1];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["{x, y} = PdfFromWeightedTrace[wtrace, weights, x];:"//OutputForm];
{x, y} = PdfFromWeightedTrace[wtrace, weights, x];
ListLinePlot[{Transpose[{x, y}]}]
Assert[AllTrue[y,#>=-10^-14 &], "PdfFromWeightedTrace returned a wrong pdf!"];
Print["========================================"]
Print["Testing BuTools function MarginalMomentsFromWeightedTrace"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
wtrace = {0.12, 1.23, 0.546, 0.6765, 1.34, 2.34};
Print["wtrace = "//OutputForm];
Print[wtrace];
weights = {12., 1., 34., 23., 8., 2.};
Print["weights = "//OutputForm];
Print[weights];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["moms = MarginalMomentsFromWeightedTrace[wtrace, weights, 3];:"//OutputForm];
moms = MarginalMomentsFromWeightedTrace[wtrace, weights, 3];
Print["moms = "//OutputForm];
Print[moms];
