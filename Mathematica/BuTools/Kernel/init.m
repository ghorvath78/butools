(* ::Package:: *)

Print["Butools V2.0"]
BuTools`Packages={"Moments","MC","RepTrans","PH"};
WriteString["stdout","Packages loaded: "];
Do[
	Get[StringJoin["BuTools`",p,"`"]];
	WriteString["stdout",p];
	WriteString["stdout","\t"];
,{p,BuTools`Packages}];
WriteString["stdout","\n"];
BuTools`Verbose=False;
BuTools`CheckInput=True;
BuTools`CheckPrecision=N[10^-12];
Print["Global variables: \nBuTools`Verbose = ",BuTools`Verbose, "\nBuTools`CheckInput = ",BuTools`CheckInput,"\nBuTools`CheckPrecision = ",BuTools`CheckPrecision];









