(* ::Package:: *)

(*
   BuTools Trace Package
*)

BeginPackage["BuTools`Trace`"];
CdfFromTrace::usage = "";
CdfFromWeightedTrace::usage = "";
IATimesFromCummulative::usage = "";
LagCorrelationsFromTrace::usage = "";
LagkJointMomentsFromTrace::usage = "";
MarginalMomentsFromTrace::usage = "";
MarginalMomentsFromWeightedTrace::usage = "";
PdfFromTrace::usage = "";
PdfFromWeightedTrace::usage = "";


Begin["`Private`"];


If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckInput"]],BuTools`CheckInput=True];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckPrecision"]],BuTools`CheckPrecision=N[10^-12]];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`Verbose"]],BuTools`Verbose=False];


CdfFromTrace[trace_]:={Sort[trace], Range[0, 1, 1./(Length[trace]-1)]};


CdfFromWeightedTrace[trace_,weights_]:=Module[{ord},
	ord=Ordering[trace];
	Return[{trace[[ord]],Accumulate[weights[[ord]]]/Total[weights]}];
];


IATimesFromCummulative[trace_]:=Differences[trace];


LagCorrelationsFromTrace[trace_, K_:3]:=
Module[{m,v,acf},
    m = Mean[trace];
    v = Variance[trace];
    Return[Table[((trace[[1;;-i-1]].trace[[i+1;;]]) / (Length[trace]-i) - m^2) / v,{i,K}]];
];


LagkJointMomentsFromTrace[trace_, K_:3, L_:1]:= 
Table[(trace[[1;;-L-1]]^i).(trace[[L+1;;]]^j)/(Length[trace]-L),{i,0,K},{j,0,K}];


MarginalMomentsFromTrace[trace_, K_:5]:=
Table[Total[trace^i],{i,1,K}]/Length[trace];


MarginalMomentsFromWeightedTrace[trace_,weights_,K_:5]:=
Table[Total[(trace^i).weights],{i,1,K}]/Total[weights];    


PdfFromTrace[trace_, intBounds_]:=
Module[{intlens,x,y,ix},
    intlens = intBounds[[2;;]] - intBounds[[1;;-2]];
    x = (intBounds[[2;;]] + intBounds[[1;;-2]]) / 2;
	ix = Range[Length[trace]];
    y = Table[Length[Select[trace,And[#>=intBounds[[i]],#<intBounds[[i+1]]]&]],{i,Length[x]}];
    y = y / intlens / Length[trace];
	Return[{x,y}];
];


PdfFromWeightedTrace[trace_, weights_, intBounds_]:=
Module[{intlens,x,y,ix},
    intlens = intBounds[[2;;]] - intBounds[[1;;-2]];
    x = (intBounds[[2;;]] + intBounds[[1;;-2]]) / 2;
	ix = Range[Length[trace]];
    y = Table[Total[weights[[Select[ix,And[trace[[#]]>=intBounds[[i]],trace[[#]]<intBounds[[i+1]]]&]]]],{i,Length[x]}];
    y = y / intlens / Total[weights];
	Return[{x,y}];
];


End[(* Private *)];
EndPackage[];
