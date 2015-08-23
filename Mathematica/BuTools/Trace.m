(* ::Package:: *)

(*
   BuTools Trace Package
*)

BeginPackage["BuTools`Trace`"];
CdfFromTrace::usage = "{x, y} = CdfFromTrace[trace]: Returns the empirical distribution function of the trace";
CdfFromWeightedTrace::usage = "{x, y} = CdfFromWeightedTrace[trace, weights]: Returns the empirical distribution function of a trace consisting of weighted data";
IATimesFromCummulative::usage = "iat = IATimesFromCummulative[trace]: Returns inter-arrival times from cummulative a trace.";
LagCorrelationsFromTrace::usage = "acf = LagCorrelationsFromTrace[trace, K]: Returns the lag-k autocorrelation of a trace";
LagkJointMomentsFromTrace::usage = "Nm = LagkJointMomentsFromTrace[trace, K, L]: Returns the lag-k joint moments of a trace";
MarginalMomentsFromTrace::usage = "moms = MarginalMomentsFromTrace[trace, K]: Returns the marginal moments of a trace";
MarginalMomentsFromWeightedTrace::usage = "moms = MarginalMomentsFromWeightedTrace[trace, weights, K]: Returns the marginal moments of a trace consisting of weighted data";
PdfFromTrace::usage = "{x, y} = PdfFromTrace[trace, intBounds]: Returns the empirical density function of a trace";
PdfFromWeightedTrace::usage = "{x, y} = PdfFromWeightedTrace[trace, weights, intBounds]: Returns the empirical density function of a trace consisting of weighted data";


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
Table[(trace^i).weights,{i,1,K}]/Total[weights];    


PdfFromTrace[trace_, intBounds_]:=
Module[{intlens,x,y,ix,str,l},
    intlens = intBounds[[2;;]] - intBounds[[1;;-2]];
    x = (intBounds[[2;;]] + intBounds[[1;;-2]]) / 2;
    str =Sort[trace]; 
	l=LengthWhile[str,#<intBounds[[1]]&];
	str=Drop[str,l];
	y={};
	Do[
		l=LengthWhile[str,#<intBounds[[i]]&];
		AppendTo[y,l];
		str=Drop[str,l];
	,{i,2,Length[intBounds]}];
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
