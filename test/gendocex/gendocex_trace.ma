ClearAll["Global`*"]
AppendTo[$Path,"/home/gabor/github/butools/Mathematica"];
<<BuTools`
Print["---BuTools: Trace package test file---"//OutputForm];
Print["Enable the verbose messages with the BuToolsVerbose flag"//OutputForm];
BuTools`Verbose = True;
Print["Enable input parameter checking with the BuToolsCheckInput flag"//OutputForm];
BuTools`CheckInput = true;
On[Assert];
tmpOut = $Output;
stream = OpenWrite["/home/gabor/github/butools/test/docex/Trace_mathematica.docex", FormatType -> InputForm, PageWidth -> Infinity];
$Output = {stream};
Unprotect[Print];
Print[args___] := Block[{$inMsg = True, result, str},
   If[MatrixQ[args],
       str = "{";
       Do[str = StringJoin[str, ToString[args[[r]], FormatType -> InputForm]]; 
            If[r < Length[args], str = StringJoin[str, ",\n "]], {r, Length[args]}];
            str = StringJoin[str, "}"];
            Print[str//OutputForm],
            result = Print[args],
            result = Print[args]
        ]] /; ! TrueQ[$inMsg];
Print["=== CdfFromTrace ==="]
Print[">>> D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};"//OutputForm];
D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};
Print[">>> D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};"//OutputForm];
D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};
Print[">>> tr = SamplesFromMAP[D0, D1, 1000000];"//OutputForm];
tr = SamplesFromMAP[D0, D1, 1000000];
Print[">>> {x, y} = CdfFromTrace[tr];"//OutputForm];
{x, y} = CdfFromTrace[tr];
Print[">>> ListLinePlot[{Transpose[{x, y}]}]"//OutputForm];
meandiff = Abs[Dot[Differences[x], 1.-y[[1;;-2]]]-Mean[tr]]/Mean[tr];
Print["=== PdfFromTrace ==="]
Print[">>> D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};"//OutputForm];
D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};
Print[">>> D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};"//OutputForm];
D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};
Print[">>> x = Range[0.0,0.5,0.01];"//OutputForm];
x = Range[0.0,0.5,0.01];
Print[">>> tr = SamplesFromMAP[D0, D1, 1000000];"//OutputForm];
tr = SamplesFromMAP[D0, D1, 1000000];
Print[">>> {x, y} = PdfFromTrace[tr, x];"//OutputForm];
{x, y} = PdfFromTrace[tr, x];
Print[">>> {a, A} = MarginalDistributionFromMAP[D0, D1];"//OutputForm];
{a, A} = MarginalDistributionFromMAP[D0, D1];
Print[">>> {xm, ym} = IntervalPdfFromPH[a, A, x];"//OutputForm];
{xm, ym} = IntervalPdfFromPH[a, A, x];
Print[">>> ListLinePlot[{Transpose[{x, y}],Transpose[{xm, ym}]}]"//OutputForm];
meandiff = Abs[Dot[x[[1;;-2]], (Differences[x]*y[[1;;-2]])]-Mean[tr]]/Mean[tr];
Print["=== MarginalMomentsFromTrace ==="]
Print[">>> D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};"//OutputForm];
D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};
Print[">>> D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};"//OutputForm];
D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};
Print[">>> tr = SamplesFromMAP[D0, D1, 1000000];"//OutputForm];
tr = SamplesFromMAP[D0, D1, 1000000];
Print[">>> moms = MarginalMomentsFromTrace[tr, 3];"//OutputForm];
moms = MarginalMomentsFromTrace[tr, 3];
Print[">>> Print[moms];"//OutputForm];
Print[moms];
Print[">>> mmoms = MarginalMomentsFromMAP[D0, D1, 3];"//OutputForm];
mmoms = MarginalMomentsFromMAP[D0, D1, 3];
Print[">>> Print[mmoms];"//OutputForm];
Print[mmoms];
Print["=== LagCorrelationsFromTrace ==="]
Print[">>> D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};"//OutputForm];
D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};
Print[">>> D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};"//OutputForm];
D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};
Print[">>> tr = SamplesFromMAP[D0, D1, 1000000];"//OutputForm];
tr = SamplesFromMAP[D0, D1, 1000000];
Print[">>> acf = LagCorrelationsFromTrace[tr, 10];"//OutputForm];
acf = LagCorrelationsFromTrace[tr, 10];
Print[">>> Print[acf];"//OutputForm];
Print[acf];
Print[">>> macf = LagCorrelationsFromMAP[D0, D1, 10];"//OutputForm];
macf = LagCorrelationsFromMAP[D0, D1, 10];
Print[">>> Print[macf];"//OutputForm];
Print[macf];
Print["=== LagkJointMomentsFromTrace ==="]
Print[">>> D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};"//OutputForm];
D0 = {{-18., 1., 4.},{2., -18., 7.},{1., 3., -32.}};
Print[">>> D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};"//OutputForm];
D1 = {{12., 1., 0.},{1., 8., 0.},{2., 1., 25.}};
Print[">>> tr = SamplesFromMAP[D0, D1, 1000000];"//OutputForm];
tr = SamplesFromMAP[D0, D1, 1000000];
Print[">>> Nm1 = LagkJointMomentsFromTrace[tr, 3, 1];"//OutputForm];
Nm1 = LagkJointMomentsFromTrace[tr, 3, 1];
Print[">>> Print[Nm1];"//OutputForm];
Print[Nm1];
Print[">>> mNm1 = LagkJointMomentsFromMAP[D0, D1, 3, 1];"//OutputForm];
mNm1 = LagkJointMomentsFromMAP[D0, D1, 3, 1];
Print[">>> Print[mNm1];"//OutputForm];
Print[mNm1];
Print[">>> Nm2 = LagkJointMomentsFromTrace[tr, 3, 2];"//OutputForm];
Nm2 = LagkJointMomentsFromTrace[tr, 3, 2];
Print[">>> Print[Nm2];"//OutputForm];
Print[Nm2];
Print[">>> mNm2 = LagkJointMomentsFromMAP[D0, D1, 3, 2];"//OutputForm];
mNm2 = LagkJointMomentsFromMAP[D0, D1, 3, 2];
Print[">>> Print[mNm2];"//OutputForm];
Print[mNm2];
Print["=== CdfFromWeightedTrace ==="]
Print[">>> wtrace = {0.12, 1.23, 0.546, 0.6765, 1.34, 2.34};"//OutputForm];
wtrace = {0.12, 1.23, 0.546, 0.6765, 1.34, 2.34};
Print[">>> weights = {12., 1., 34., 23., 8., 2.};"//OutputForm];
weights = {12., 1., 34., 23., 8., 2.};
Print[">>> {x, y} = CdfFromWeightedTrace[wtrace, weights];"//OutputForm];
{x, y} = CdfFromWeightedTrace[wtrace, weights];
Print[">>> ListLinePlot[{Transpose[{x, y}]}]"//OutputForm];
Print["=== PdfFromWeightedTrace ==="]
Print[">>> wtrace = {0.12, 1.23, 0.546, 0.6765, 1.34, 2.34};"//OutputForm];
wtrace = {0.12, 1.23, 0.546, 0.6765, 1.34, 2.34};
Print[">>> weights = {12., 1., 34., 23., 8., 2.};"//OutputForm];
weights = {12., 1., 34., 23., 8., 2.};
Print[">>> x = Range[0.0,3.0,0.1];"//OutputForm];
x = Range[0.0,3.0,0.1];
Print[">>> {x, y} = PdfFromWeightedTrace[wtrace, weights, x];"//OutputForm];
{x, y} = PdfFromWeightedTrace[wtrace, weights, x];
Print[">>> ListLinePlot[{Transpose[{x, y}]}]"//OutputForm];
Print["=== MarginalMomentsFromWeightedTrace ==="]
Print[">>> wtrace = {0.12, 1.23, 0.546, 0.6765, 1.34, 2.34};"//OutputForm];
wtrace = {0.12, 1.23, 0.546, 0.6765, 1.34, 2.34};
Print[">>> weights = {12., 1., 34., 23., 8., 2.};"//OutputForm];
weights = {12., 1., 34., 23., 8., 2.};
Print[">>> moms = MarginalMomentsFromWeightedTrace[wtrace, weights, 3];"//OutputForm];
moms = MarginalMomentsFromWeightedTrace[wtrace, weights, 3];
Print[">>> Print[moms];"//OutputForm];
Print[moms];
$Output = tmpOut;Close[stream];

