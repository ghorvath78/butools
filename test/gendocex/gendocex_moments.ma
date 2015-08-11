ClearAll["Global`*"]
AppendTo[$Path,"/home/gabor/github/butools/Mathematica"];
<<BuTools`
Print["---BuTools: Moments package test file---"//OutputForm];
Print["Enable the verbose messages with the BuToolsVerbose flag"//OutputForm];
BuTools`Verbose = True;
Print["Enable input parameter checking with the BuToolsCheckInput flag"//OutputForm];
BuTools`CheckInput = true;
On[Assert];
tmpOut = $Output;
stream = OpenWrite["/home/gabor/github/butools/test/docex/Moments_mathematica.docex", FormatType -> InputForm, PageWidth -> Infinity];
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
Print["=== NormMomsFromMoms ==="]
Print[">>> M = {1.2, 5., 38., 495., 9215.};"//OutputForm];
M = {1.2, 5., 38., 495., 9215.};
Print[">>> nmoms = NormMomsFromMoms[M];"//OutputForm];
nmoms = NormMomsFromMoms[M];
Print[">>> Print[nmoms];"//OutputForm];
Print[nmoms];
Print[">>> moms = MomsFromNormMoms[nmoms];"//OutputForm];
moms = MomsFromNormMoms[nmoms];
Print[">>> Print[moms];"//OutputForm];
Print[moms];
Print[">>> err = Norm[moms-M];"//OutputForm];
err = Norm[moms-M];
Print[">>> Print[err];"//OutputForm];
Print[err];
Print["=== MomsFromNormMoms ==="]
Print[">>> M = {1.2, 5., 38., 495., 9215.};"//OutputForm];
M = {1.2, 5., 38., 495., 9215.};
Print[">>> nmoms = NormMomsFromMoms[M];"//OutputForm];
nmoms = NormMomsFromMoms[M];
Print[">>> Print[nmoms];"//OutputForm];
Print[nmoms];
Print[">>> moms = MomsFromNormMoms[nmoms];"//OutputForm];
moms = MomsFromNormMoms[nmoms];
Print[">>> Print[moms];"//OutputForm];
Print[moms];
Print[">>> err = Norm[moms-M];"//OutputForm];
err = Norm[moms-M];
Print[">>> Print[err];"//OutputForm];
Print[err];
Print["=== ReducedMomsFromMoms ==="]
Print[">>> M = {1.2, 5., 38., 495., 9215.};"//OutputForm];
M = {1.2, 5., 38., 495., 9215.};
Print[">>> rmoms = ReducedMomsFromMoms[M];"//OutputForm];
rmoms = ReducedMomsFromMoms[M];
Print[">>> Print[rmoms];"//OutputForm];
Print[rmoms];
Print[">>> moms = MomsFromReducedMoms[rmoms];"//OutputForm];
moms = MomsFromReducedMoms[rmoms];
Print[">>> Print[moms];"//OutputForm];
Print[moms];
Print[">>> err = Norm[moms-M];"//OutputForm];
err = Norm[moms-M];
Print[">>> Print[err];"//OutputForm];
Print[err];
Print["=== MomsFromReducedMoms ==="]
Print[">>> M = {1.2, 5., 38., 495., 9215.};"//OutputForm];
M = {1.2, 5., 38., 495., 9215.};
Print[">>> rmoms = ReducedMomsFromMoms[M];"//OutputForm];
rmoms = ReducedMomsFromMoms[M];
Print[">>> Print[rmoms];"//OutputForm];
Print[rmoms];
Print[">>> moms = MomsFromReducedMoms[rmoms];"//OutputForm];
moms = MomsFromReducedMoms[rmoms];
Print[">>> Print[moms];"//OutputForm];
Print[moms];
Print[">>> err = Norm[moms-M];"//OutputForm];
err = Norm[moms-M];
Print[">>> Print[err];"//OutputForm];
Print[err];
Print["=== FactorialMomsFromMoms ==="]
Print[">>> M = {1.3, 2.4, 6.03, 20.5, 89.5, 474.9};"//OutputForm];
M = {1.3, 2.4, 6.03, 20.5, 89.5, 474.9};
Print[">>> fmoms = FactorialMomsFromMoms[M];"//OutputForm];
fmoms = FactorialMomsFromMoms[M];
Print[">>> Print[fmoms];"//OutputForm];
Print[fmoms];
Print[">>> moms = MomsFromFactorialMoms[fmoms];"//OutputForm];
moms = MomsFromFactorialMoms[fmoms];
Print[">>> Print[moms];"//OutputForm];
Print[moms];
Print[">>> err = Norm[moms-M];"//OutputForm];
err = Norm[moms-M];
Print[">>> Print[err];"//OutputForm];
Print[err];
Print["=== MomsFromFactorialMoms ==="]
Print[">>> M = {1.3, 2.4, 6.03, 20.5, 89.5, 474.9};"//OutputForm];
M = {1.3, 2.4, 6.03, 20.5, 89.5, 474.9};
Print[">>> fmoms = FactorialMomsFromMoms[M];"//OutputForm];
fmoms = FactorialMomsFromMoms[M];
Print[">>> Print[fmoms];"//OutputForm];
Print[fmoms];
Print[">>> moms = MomsFromFactorialMoms[fmoms];"//OutputForm];
moms = MomsFromFactorialMoms[fmoms];
Print[">>> Print[moms];"//OutputForm];
Print[moms];
Print[">>> err = Norm[moms-M];"//OutputForm];
err = Norm[moms-M];
Print[">>> Print[err];"//OutputForm];
Print[err];
Print["=== JFactorialMomsFromJMoms ==="]
Print[">>> MM = {{0.7, 2., 3., 4.},{5., 6., 7., 8.},{9., 10., 11., 12.}};"//OutputForm];
MM = {{0.7, 2., 3., 4.},{5., 6., 7., 8.},{9., 10., 11., 12.}};
Print[">>> JFmoms = JFactorialMomsFromJMoms[MM];"//OutputForm];
JFmoms = JFactorialMomsFromJMoms[MM];
Print[">>> Print[JFmoms];"//OutputForm];
Print[JFmoms];
Print[">>> Jmoms = JMomsFromJFactorialMoms[JFmoms];"//OutputForm];
Jmoms = JMomsFromJFactorialMoms[JFmoms];
Print[">>> Print[Jmoms];"//OutputForm];
Print[Jmoms];
Print[">>> err = Norm[Jmoms-MM];"//OutputForm];
err = Norm[Jmoms-MM];
Print[">>> Print[err];"//OutputForm];
Print[err];
Print["=== JMomsFromJFactorialMoms ==="]
Print[">>> MM = {{0.7, 2., 3., 4.},{5., 6., 7., 8.},{9., 10., 11., 12.}};"//OutputForm];
MM = {{0.7, 2., 3., 4.},{5., 6., 7., 8.},{9., 10., 11., 12.}};
Print[">>> JFmoms = JFactorialMomsFromJMoms[MM];"//OutputForm];
JFmoms = JFactorialMomsFromJMoms[MM];
Print[">>> Print[JFmoms];"//OutputForm];
Print[JFmoms];
Print[">>> Jmoms = JMomsFromJFactorialMoms[JFmoms];"//OutputForm];
Jmoms = JMomsFromJFactorialMoms[JFmoms];
Print[">>> Print[Jmoms];"//OutputForm];
Print[Jmoms];
Print[">>> err = Norm[Jmoms-MM];"//OutputForm];
err = Norm[Jmoms-MM];
Print[">>> Print[err];"//OutputForm];
Print[err];
Print["=== HankelMomsFromMoms ==="]
Print[">>> M = {1.3, 2.4, 6.03, 20.5, 89.5, 474.9};"//OutputForm];
M = {1.3, 2.4, 6.03, 20.5, 89.5, 474.9};
Print[">>> hmoms = HankelMomsFromMoms[M];"//OutputForm];
hmoms = HankelMomsFromMoms[M];
Print[">>> Print[hmoms];"//OutputForm];
Print[hmoms];
Print[">>> moms = MomsFromHankelMoms[hmoms];"//OutputForm];
moms = MomsFromHankelMoms[hmoms];
Print[">>> Print[moms];"//OutputForm];
Print[moms];
Print[">>> err = Norm[moms-M];"//OutputForm];
err = Norm[moms-M];
Print[">>> Print[err];"//OutputForm];
Print[err];
Print["=== MomsFromHankelMoms ==="]
Print[">>> M = {1.3, 2.4, 6.03, 20.5, 89.5, 474.9};"//OutputForm];
M = {1.3, 2.4, 6.03, 20.5, 89.5, 474.9};
Print[">>> hmoms = HankelMomsFromMoms[M];"//OutputForm];
hmoms = HankelMomsFromMoms[M];
Print[">>> Print[hmoms];"//OutputForm];
Print[hmoms];
Print[">>> moms = MomsFromHankelMoms[hmoms];"//OutputForm];
moms = MomsFromHankelMoms[hmoms];
Print[">>> Print[moms];"//OutputForm];
Print[moms];
Print[">>> err = Norm[moms-M];"//OutputForm];
err = Norm[moms-M];
Print[">>> Print[err];"//OutputForm];
Print[err];
Print["=== CheckMoments ==="]
Print[">>> M = {1.2, 5., 8., 29., 3412.};"//OutputForm];
M = {1.2, 5., 8., 29., 3412.};
Print[">>> flag = CheckMoments[M];"//OutputForm];
flag = CheckMoments[M];
Print[">>> Print[flag];"//OutputForm];
Print[flag];
Print[">>> M = {1.3, 2.4, 6.03, 20.5, 89.5};"//OutputForm];
M = {1.3, 2.4, 6.03, 20.5, 89.5};
Print[">>> flag = CheckMoments[M];"//OutputForm];
flag = CheckMoments[M];
Print[">>> Print[flag];"//OutputForm];
Print[flag];
$Output = tmpOut;Close[stream];

