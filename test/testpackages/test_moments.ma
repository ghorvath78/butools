ClearAll["Global`*"]
AppendTo[$Path,"/home/gabor/github/butools/Mathematica"];
<<BuTools`
Print["---BuTools: Moments package test file---"//OutputForm];
Print["Enable the verbose messages with the BuToolsVerbose flag"//OutputForm];
BuTools`Verbose = True;
Print["Enable input parameter checking with the BuToolsCheckInput flag"//OutputForm];
BuTools`CheckInput = true;
On[Assert];
Print["========================================"]
Print["Testing BuTools function NormMomsFromMoms"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
M = {1.2, 5., 38., 495., 9215.};
Print["M = "//OutputForm];
Print[M];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["nmoms = NormMomsFromMoms[M];:"//OutputForm];
nmoms = NormMomsFromMoms[M];
Print["nmoms = "//OutputForm];
Print[nmoms];
Print["moms = MomsFromNormMoms[nmoms];:"//OutputForm];
moms = MomsFromNormMoms[nmoms];
Print["moms = "//OutputForm];
Print[moms];
Print["err = Norm[moms-M];:"//OutputForm];
err = Norm[moms-M];
Print["err = "//OutputForm];
Print[err];
Assert[err<10^-10, "Calling the moment conversion and its inverse did not give back the original moments!"];
Print["========================================"]
Print["Testing BuTools function MomsFromNormMoms"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
M = {1.2, 5., 38., 495., 9215.};
Print["M = "//OutputForm];
Print[M];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["nmoms = NormMomsFromMoms[M];:"//OutputForm];
nmoms = NormMomsFromMoms[M];
Print["nmoms = "//OutputForm];
Print[nmoms];
Print["moms = MomsFromNormMoms[nmoms];:"//OutputForm];
moms = MomsFromNormMoms[nmoms];
Print["moms = "//OutputForm];
Print[moms];
Print["err = Norm[moms-M];:"//OutputForm];
err = Norm[moms-M];
Print["err = "//OutputForm];
Print[err];
Assert[err<10^-10, "Calling the moment conversion and its inverse did not give back the original moments!"];
Print["========================================"]
Print["Testing BuTools function ReducedMomsFromMoms"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
M = {1.2, 5., 38., 495., 9215.};
Print["M = "//OutputForm];
Print[M];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["rmoms = ReducedMomsFromMoms[M];:"//OutputForm];
rmoms = ReducedMomsFromMoms[M];
Print["rmoms = "//OutputForm];
Print[rmoms];
Print["moms = MomsFromReducedMoms[rmoms];:"//OutputForm];
moms = MomsFromReducedMoms[rmoms];
Print["moms = "//OutputForm];
Print[moms];
Print["err = Norm[moms-M];:"//OutputForm];
err = Norm[moms-M];
Print["err = "//OutputForm];
Print[err];
Assert[err<10^-10, "Calling the moment conversion and its inverse did not give back the original moments!"];
Print["========================================"]
Print["Testing BuTools function MomsFromReducedMoms"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
M = {1.2, 5., 38., 495., 9215.};
Print["M = "//OutputForm];
Print[M];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["rmoms = ReducedMomsFromMoms[M];:"//OutputForm];
rmoms = ReducedMomsFromMoms[M];
Print["rmoms = "//OutputForm];
Print[rmoms];
Print["moms = MomsFromReducedMoms[rmoms];:"//OutputForm];
moms = MomsFromReducedMoms[rmoms];
Print["moms = "//OutputForm];
Print[moms];
Print["err = Norm[moms-M];:"//OutputForm];
err = Norm[moms-M];
Print["err = "//OutputForm];
Print[err];
Assert[err<10^-10, "Calling the moment conversion and its inverse did not give back the original moments!"];
Print["========================================"]
Print["Testing BuTools function FactorialMomsFromMoms"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
M = {1.3, 2.4, 6.03, 20.5, 89.5, 474.9};
Print["M = "//OutputForm];
Print[M];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["fmoms = FactorialMomsFromMoms[M];:"//OutputForm];
fmoms = FactorialMomsFromMoms[M];
Print["fmoms = "//OutputForm];
Print[fmoms];
Print["moms = MomsFromFactorialMoms[fmoms];:"//OutputForm];
moms = MomsFromFactorialMoms[fmoms];
Print["moms = "//OutputForm];
Print[moms];
Print["err = Norm[moms-M];:"//OutputForm];
err = Norm[moms-M];
Print["err = "//OutputForm];
Print[err];
Assert[err<10^-10, "Calling the moment conversion and its inverse did not give back the original moments!"];
Print["========================================"]
Print["Testing BuTools function MomsFromFactorialMoms"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
M = {1.3, 2.4, 6.03, 20.5, 89.5, 474.9};
Print["M = "//OutputForm];
Print[M];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["fmoms = FactorialMomsFromMoms[M];:"//OutputForm];
fmoms = FactorialMomsFromMoms[M];
Print["fmoms = "//OutputForm];
Print[fmoms];
Print["moms = MomsFromFactorialMoms[fmoms];:"//OutputForm];
moms = MomsFromFactorialMoms[fmoms];
Print["moms = "//OutputForm];
Print[moms];
Print["err = Norm[moms-M];:"//OutputForm];
err = Norm[moms-M];
Print["err = "//OutputForm];
Print[err];
Assert[err<10^-10, "Calling the moment conversion and its inverse did not give back the original moments!"];
Print["========================================"]
Print["Testing BuTools function JFactorialMomsFromJMoms"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
MM = {{0.7, 2., 3., 4.},{5., 6., 7., 8.},{9., 10., 11., 12.}};
Print["MM = "//OutputForm];
Print[MM];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["JFmoms = JFactorialMomsFromJMoms[MM];:"//OutputForm];
JFmoms = JFactorialMomsFromJMoms[MM];
Print["JFmoms = "//OutputForm];
Print[JFmoms];
Print["Jmoms = JMomsFromJFactorialMoms[JFmoms];:"//OutputForm];
Jmoms = JMomsFromJFactorialMoms[JFmoms];
Print["Jmoms = "//OutputForm];
Print[Jmoms];
Print["err = Norm[Jmoms-MM];:"//OutputForm];
err = Norm[Jmoms-MM];
Print["err = "//OutputForm];
Print[err];
Assert[err<10^-10, "Calling the moment conversion and its inverse did not give back the original moments!"];
Print["========================================"]
Print["Testing BuTools function JMomsFromJFactorialMoms"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
MM = {{0.7, 2., 3., 4.},{5., 6., 7., 8.},{9., 10., 11., 12.}};
Print["MM = "//OutputForm];
Print[MM];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["JFmoms = JFactorialMomsFromJMoms[MM];:"//OutputForm];
JFmoms = JFactorialMomsFromJMoms[MM];
Print["JFmoms = "//OutputForm];
Print[JFmoms];
Print["Jmoms = JMomsFromJFactorialMoms[JFmoms];:"//OutputForm];
Jmoms = JMomsFromJFactorialMoms[JFmoms];
Print["Jmoms = "//OutputForm];
Print[Jmoms];
Print["err = Norm[Jmoms-MM];:"//OutputForm];
err = Norm[Jmoms-MM];
Print["err = "//OutputForm];
Print[err];
Assert[err<10^-10, "Calling the moment conversion and its inverse did not give back the original moments!"];
Print["========================================"]
Print["Testing BuTools function HankelMomsFromMoms"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
M = {1.3, 2.4, 6.03, 20.5, 89.5, 474.9};
Print["M = "//OutputForm];
Print[M];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["hmoms = HankelMomsFromMoms[M];:"//OutputForm];
hmoms = HankelMomsFromMoms[M];
Print["hmoms = "//OutputForm];
Print[hmoms];
Print["moms = MomsFromHankelMoms[hmoms];:"//OutputForm];
moms = MomsFromHankelMoms[hmoms];
Print["moms = "//OutputForm];
Print[moms];
Print["err = Norm[moms-M];:"//OutputForm];
err = Norm[moms-M];
Print["err = "//OutputForm];
Print[err];
Assert[err<10^-10, "Calling the moment conversion and its inverse did not give back the original moments!"];
Print["========================================"]
Print["Testing BuTools function MomsFromHankelMoms"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
M = {1.3, 2.4, 6.03, 20.5, 89.5, 474.9};
Print["M = "//OutputForm];
Print[M];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["hmoms = HankelMomsFromMoms[M];:"//OutputForm];
hmoms = HankelMomsFromMoms[M];
Print["hmoms = "//OutputForm];
Print[hmoms];
Print["moms = MomsFromHankelMoms[hmoms];:"//OutputForm];
moms = MomsFromHankelMoms[hmoms];
Print["moms = "//OutputForm];
Print[moms];
Print["err = Norm[moms-M];:"//OutputForm];
err = Norm[moms-M];
Print["err = "//OutputForm];
Print[err];
Assert[err<10^-10, "Calling the moment conversion and its inverse did not give back the original moments!"];
Print["========================================"]
Print["Testing BuTools function CheckMoments"]
Print["Input:"//OutputForm];
Print["------"//OutputForm];
M = {1.2, 5., 8., 29., 3412.};
Print["M = "//OutputForm];
Print[M];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["flag = CheckMoments[M];:"//OutputForm];
flag = CheckMoments[M];
Print["flag = "//OutputForm];
Print[flag];
Assert[flag==False, "CheckMoments did not recognize a valid moment sequence!"];
Print["Input:"//OutputForm];
Print["------"//OutputForm];
M = {1.3, 2.4, 6.03, 20.5, 89.5};
Print["M = "//OutputForm];
Print[M];
Print["Test:"//OutputForm];
Print["-----"//OutputForm];
Print["flag = CheckMoments[M];:"//OutputForm];
flag = CheckMoments[M];
Print["flag = "//OutputForm];
Print[flag];
Assert[flag==True, "CheckMoments did not recognize an invalid moment sequence!"];