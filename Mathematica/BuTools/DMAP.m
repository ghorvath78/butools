(* ::Package:: *)

(*
   BuTools DMAP Package
*)


BeginPackage["BuTools`DMAP`"];
CanonicalFromDMAP2::usage = "{G0, G1} = CanonicalFromDMAP2[D0, D1, prec]: Returns the canonical form of an order-2 discrete Markovian arrival process.";
CheckDMAPRepresentation::usage = "r = CheckDMAPRepresentation[D0, D1, prec]: Checks if the input matrixes define a discrete time MAP.";
CheckDMMAPRepresentation::usage = "r = CheckDMMAPRepresentation[D, prec]: Checks if the input matrixes define a discrete time MMAP.";
CheckDRAPRepresentation::usage = "r = CheckDRAPRepresentation[H, prec]: Checks if the input matrixes define a discrete time RAP.";
CheckDMRAPRepresentation::usage = "r = CheckDMRAPRepresentation[H, prec]: Checks if the input matrixes define a discrete time MRAP.";
LagCorrelationsFromDRAP::usage = "acf = LagCorrelationsFromDRAP[H0, H1, L, prec]: Returns the lag autocorrelations of a discrete rational arrival process.";
LagCorrelationsFromDMAP::usage = "acf = LagCorrelationsFromDMAP[D0, D1, L, prec]: Returns the lag autocorrelations of a discrete Markovian arrival process.";
LagkJointMomentsFromDRAP::usage = "Nm = LagkJointMomentsFromDRAP[H0, H1, K, L, prec]: Returns the lag-k joint moments of a discrete rational arrival process.";
LagkJointMomentsFromDMRAP::usage = "Nm = LagkJointMomentsFromDMRAP[H, K, L, prec]: Returns the lag-k joint moments of a discrete marked rational arrival process.";
LagkJointMomentsFromDMAP::usage = "Nm = LagkJointMomentsFromDMAP[D0, D1, K, L, prec]: Returns the lag-k joint moments of a discrete Markovian arrival process.";
LagkJointMomentsFromDMMAP::usage = "Nm = LagkJointMomentsFromDMMAP[D, K, L, prec]: Returns the lag-k joint moments of a discrete marked Markovian arrival process.";
DMAP2FromMoments::usage = "{D0, D1} = DMAP2FromMoments[moms, corr1]: Returns a DMAP(2) which has the same 3 marginal moments and lag-1 autocorrelation as given.";
DMAPFromDRAP::usage = "{D0, D1} = DMAPFromDRAP[H0, H1, precision]: Obtains a Markovian representation of a discrete rational arrival process of the same size, if possible.";
DMMAPFromDMRAP::usage = "D = DMMAPFromDMRAP[H, precision]: Obtains a Markovian representation of a discrete marked rational arrival process of the same size, if possible.";
MarginalDistributionFromDRAP::usage = "{alpha, A} = MarginalDistributionFromDRAP[H0, H1, precision]: Returns the matrix-geometrically distributed marginal of a discrete rational arrival process.";
MarginalDistributionFromDMRAP::usage = "{alpha, A} = MarginalDistributionFromDMRAP[H, precision]: Returns the matrix-geometrically distributed marginal of a discrete marked rational arrival process.";
MarginalDistributionFromDMAP::usage = "{alpha, A} = MarginalDistributionFromDMAP[D0, D1, precision]: Returns the discrete phase type distributed marginal distribution of a discrete Markovian arrival process.";
MarginalDistributionFromDMMAP::usage = "{alpha, A} = MarginalDistributionFromDMMAP[D, precision]: Returns the discrete phase type distributed marginal of a discrete marked Markovian arrival process.     ";
MarginalMomentsFromDRAP::usage = "moms = MarginalMomentsFromDRAP[H0, H1, K, precision]: Returns the moments of the marginal distribution of a discrete rational arrival process.";
MarginalMomentsFromDMRAP::usage = "moms = MarginalMomentsFromDMRAP[H, K, precision]: Returns the moments of the marginal distribution of a discrete marked rational arrival process.";
MarginalMomentsFromDMAP::usage = "moms = MarginalMomentsFromDMAP[D0, D1, K, precision]: Returns the moments of the marginal distribution of a discrete Markovian arrival process.";
MarginalMomentsFromDMMAP::usage = "moms = MarginalMomentsFromDMMAP[D, K, precision]: Returns the moments of the marginal distribution of a discrete marked Markovian arrival process.     ";
DRAPFromMoments::usage = "{H0, H1} = DRAPFromMoments[moms, Nm]: Creates a discrete rational arrival process that has the same marginal and lag-1 joint moments as given.";
DMRAPFromMoments::usage = "H = DMRAPFromMoments[moms, Nm]: Creates a discrete marked rational arrival process that has the same marginal and lag-1 joint moments as given.";
RandomDMAP::usage = "{D0, D1} = RandomDMAP[order, mean, zeroEntries, maxTrials, prec]: Returns a random discrete Markovian arrival process.";
RandomDMMAP::usage = "D = RandomDMMAP[order, types, mean, zeroEntries, maxTrials, prec]: Returns a random discrete marked Markovian arrival process.";
SamplesFromDMAP::usage = "x = SamplesFromDMAP[D0, D1, K, prec]: Generates random samples from a discrete Markovian arrival process.";
SamplesFromDMMAP::usage = "x = SamplesFromDMMAP[D, K, prec]: Generates random samples from a discrete marked Markovian arrival process.";


Begin["`Private`"];


Needs["BuTools`RepTrans`"];
Needs["BuTools`Moments`"];
Needs["BuTools`MC`"];
Needs["BuTools`DPH`"];
Needs["BuTools`MAP`"];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckInput"]],BuTools`CheckInput=True];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckPrecision"]],BuTools`CheckPrecision=N[10^-12]];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`Verbose"]],BuTools`Verbose=False];


CanonicalFromDMAP2[D0_,D1_]:=
Module[ {s1,s2,\[Gamma], g0,g1,\[Alpha]v,W,w1,w2,a1,a,b},
If[BuTools`CheckInput && Not[Dimensions[D0][[1]] ==2],
If[BuTools`Verbose,Throw["DMAP2Canonical: Size is not 2!"]];
];

If[BuTools`CheckInput && !CheckDMAPRepresentation[D0,D1],Throw["DMAP2Canonical: Input is not a valid DMAP representation!"]];

{s1,s2}=Sort[Eigenvalues[D0],Abs[#1]>Abs[#2]||(Abs[#1]==Abs[#2] && Re[#1]>Re[#2])||(Abs[#1]==Abs[#2]&&Re[#1]==Re[#2]&&Im[#1]>=Im[#2])&];

If[s2>=0,
{g0,g1}=CanonicalFromMAP2[D0-IdentityMatrix[2],D1];
Return[{g0+IdentityMatrix[2],g1}];
];

\[Alpha]v=DRPSolve[Inverse[IdentityMatrix[2]-D0].D1];
\[Gamma]=Eigenvalues[Inverse[IdentityMatrix[2]-D0].D1];
\[Gamma]=\[Gamma][[2]];

w1=1/(s1-s2) (D0.{1,1}-s2{1,1});
w2={1,1}-w1;
W=Transpose[{w1,w2}];
a1=(1-s1)(\[Alpha]v.W)[[1]];

If[\[Gamma]>= 0,
  (* Positive correlation *)
{g0,g1}={\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{
RowBox[{"s1", "+", "s2"}], 
RowBox[{"a", " ", 
RowBox[{"(", 
RowBox[{"1", "-", "s1", "-", "s2"}], ")"}]}]},
{
FractionBox[
RowBox[{"s1", " ", "s2"}], 
RowBox[{"a", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", "s1", "+", "s2"}], ")"}]}]], "0"}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", "a"}], ")"}], 
RowBox[{"(", 
RowBox[{"1", "-", "s1", "-", "s2"}], ")"}], " "}], "0"},
{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
FractionBox[
RowBox[{"s1", " ", "s2"}], 
RowBox[{"a", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", "s1", "+", "s2"}], ")"}]}]]}], ")"}], "b"}], 
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
FractionBox[
RowBox[{"s1", " ", "s2"}], 
RowBox[{"a", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", "s1", "+", "s2"}], ")"}]}]]}], ")"}], 
RowBox[{"(", 
RowBox[{"1", "-", "b"}], ")"}]}]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\)}/.b->1+(a (-1+s1+s2-s1 s2) \[Gamma])/((-1+a) (-s1 s2+a (-1+s1+s2)))/.
a->-(1/(2 (-1+s1) (-1+s1+s2)^2))(1-4 s1+a1 s1+5 s1^2-a1 s1^2-2 s1^3-2 s2-a1 s2+5 s1 s2-3 s1^2 s2+s2^2+a1 s2^2-s1 s2^2-\[Gamma]+3 s1 \[Gamma]-a1 s1 \[Gamma]-3 s1^2 \[Gamma]+a1 s1^2 \[Gamma]+s1^3 \[Gamma]+s2 \[Gamma]+a1 s2 \[Gamma]-2 s1 s2 \[Gamma]+s1^2 s2 \[Gamma]-a1 s2^2 \[Gamma]+\[Sqrt]((-1+s1+s2)^2 ((-1+s1^2 (-2+\[Gamma])+\[Gamma]+s2 (1+a1-a1 \[Gamma])+s1 (3-a1-s2-2 \[Gamma]+a1 \[Gamma]))^2-4 (-1+s1) (-s1^3 (-1+\[Gamma])+a1 (-1+s2) s2 (-1+\[Gamma])+s1^2 (-2+a1+s2+2 \[Gamma]-a1 \[Gamma])+s1 (1-a1-s2-\[Gamma]+a1 \[Gamma])))));
Return[{g0,g1}];
];

If[\[Gamma]< 0,
  (* Negativ correlation *) 
{g0,g1}={\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{
RowBox[{"s1", "+", "s2"}], 
RowBox[{"a", " ", 
RowBox[{"(", 
RowBox[{"1", "-", "s1", "-", "s2"}], ")"}]}]},
{
FractionBox[
RowBox[{"s1", " ", "s2"}], 
RowBox[{"a", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", "s1", "+", "s2"}], ")"}]}]], "0"}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{
RowBox[{"0", " "}], 
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", "a"}], ")"}], 
RowBox[{"(", 
RowBox[{"1", "-", "s1", "-", "s2"}], ")"}]}]},
{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
FractionBox[
RowBox[{"s1", " ", "s2"}], 
RowBox[{"a", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", "s1", "+", "s2"}], ")"}]}]]}], ")"}], "b"}], 
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
FractionBox[
RowBox[{"s1", " ", "s2"}], 
RowBox[{"a", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", "s1", "+", "s2"}], ")"}]}]]}], ")"}], 
RowBox[{"(", 
RowBox[{"1", "-", "b"}], ")"}]}]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\)}/.b->-((a (-1+s1) (-1+s2) \[Gamma])/((-1+a) (-a+a s1+a s2-s1 s2)))/.
a->(a1 s1-a1 s1^2+s2-a1 s2-3 s1 s2+2 s1^2 s2-s2^2+a1 s2^2+s1 s2^2+s1 \[Gamma]-a1 s1 \[Gamma]-2 s1^2 \[Gamma]+a1 s1^2 \[Gamma]+s1^3 \[Gamma]+a1 s2 \[Gamma]-a1 s2^2 \[Gamma]+\[Sqrt](-4 (-1+s1) s1 s2 (-1+s1+s2) (a1 (s1-s2) (-1+\[Gamma])+(-1+s1) (s2+(-1+s1) \[Gamma]))+(a1 (-s1+s1^2+s2-s2^2) (-1+\[Gamma])+(-1+s1) ((-1+2 s1) s2+s2^2+(-1+s1) s1 \[Gamma]))^2))/(2 (-1+s1+s2) (a1 (s1-s2) (-1+\[Gamma])+(-1+s1) (s2+(-1+s1) \[Gamma])));
Return[{g0,g1}];
];
];


CheckDMAPRepresentation[D0_,D1_,precision_:Null] :=
Module[{prec},
	If[Not[NumericQ[precision]], prec=BuTools`CheckPrecision, prec=precision];
	If[Not[CheckProbMatrix[D0,True,prec]],Return[False]];
	If[Dimensions[D0]!=Dimensions[D1],
		If[BuTools`Verbose,Print["CheckDMAPRepresentation: D0 and D1 have different sizes!"]]; 
		Return[False]
	];
	If[ Min[D0,D1]<-prec,
		If[BuTools`Verbose,Print["CheckDMAPRepresentation: D0 and/or D1 has negative element!"]];
		Return[False]
	];
	If[Max[Abs[1-Total[D0+D1,{2}]]]>prec,
		If[BuTools`Verbose==True,Print["CheckMAPRepresentation: A rowsum of D0+D1 is not 1 (at precision ",prec,")!"]];
		Return[False]
	];
	Return[True]
];


CheckDMMAPRepresentation[D_,precision_:Null]:=
Module[{prec},
	If[Not[NumericQ[precision]], prec=BuTools`CheckPrecision, prec=precision];
	If[Min[D]<-prec,
		If[BuTools`Verbose,Print["CheckDMMAPRepresentation: One of the matrices has a negative element!"]];
		Return[False]
	];
	Return[CheckDMAPRepresentation[D[[1]],Sum[D[[i]],{i,2,Length[D]}],prec]];
];


CheckDRAPRepresentation[D0_,D1_, precision_:Null]:=
Module[ {h,size1,size2,size,eig,ceig,reig,prec},
	If[Not[NumericQ[precision]], prec=BuTools`CheckPrecision, prec=precision];
	If[ Dimensions[D0][[1]]!=Dimensions[D0][[2]],
		If[BuTools`Verbose,Print["CheckDRAPRepresentation: D0 is not a quadratic matrix!"]];
		Return[False]
	];
	If[Dimensions[D1][[1]]!=Dimensions[D1][[2]],
		If[BuTools`Verbose,Print["CheckDRAPRepresentation: D1 is not a quadratic matrix!"]];
		Return[False]
	];
	If[Dimensions[D0]!=Dimensions[D1],
		If[BuTools`Verbose,Print["CheckDRAPRepresentation: D0 and D1 have different sizes!"]]; 
		Return[False]
	];
	If[Max[Abs[1-Total[D0+D1,{2}]]] > prec,
		If[BuTools`Verbose,Print["CheckDRAPRepresentation: A rowsum of D0+D1 is not 1! (precision:",prec,")"]];
		Return[False]
	];
	
	eig=Sort[Eigenvalues[D0],Abs[#1]>Abs[#2]||(Abs[#1]==Abs[#2] && Re[#1]>Re[#2])||(Abs[#1]==Abs[#2]&&Re[#1]==Re[#2]&&Im[#1]>=Im[#2])&];

	If[Abs[Im[eig[[1]]]] > prec,
		If[BuTools`Verbose,Print["CheckDRAPRepresentation: The dominant eigenvalue of D0 is complex!"]];
		Return[False]
	];

	If[eig[[1]]>1+prec,
		If[BuTools`Verbose,Print["CheckDRAPRepresentation: The dominant eigenvalue of D0 is greater than 1!"]];
		Return[False]
	];

	If[Abs[eig[[1]]]==Abs[eig[[2]]],
		If[BuTools`Verbose,Print["CheckDRAPRepresentation Warning: D0 has more than one eigenvalue with the same absolute value as the dominant eigenvalue!"]];
	];
	Return[True]
];


CheckDMRAPRepresentation[H_, precision_:Null] := 
	CheckDRAPRepresentation[H[[1]],Sum[H[[i]],{i,2,Length[H]}],precision];


LagCorrelationsFromDRAP[H0_,H1_,L_:1]:=
Module[{H0i,P,pi,moms,acf},
	If[BuTools`CheckInput && !CheckDRAPRepresentation[H0,H1],Throw["LagCorrelationsFromDRAP: Input is not a valid DRAP representation!"]];
	H0i = Inverse[IdentityMatrix[Dimensions[H0][[1]]]-H0];
	P = H0i.H1;
	pi = DRPSolve[P];
	moms = MomentsFromMG[pi,H0,2];
	pi = pi.H0i.P;
	acf = {};
	Do[
		AppendTo[acf,(Total[pi.H0i]-moms[[1]]^2) / (moms[[2]]-moms[[1]]^2)];
		pi = pi.P;
	,{L}];
	Return[acf];
];


LagCorrelationsFromDMAP[D0_,D1_,L_:1]:=
(If[BuTools`CheckInput && !CheckDMAPRepresentation[D0,D1],Throw["LagCorrelationsFromDMAP: Input is not a valid DMAP representation!"]];
Return[LagCorrelationsFromDRAP[D0,D1,L]])


LagkJointMomentsFromDRAP[H0_,H1_,K_:0, L_:1]:=
(If[BuTools`CheckInput && !CheckDRAPRepresentation[H0,H1],Throw["LagkJointMomentsFromDRAP: Input is not a valid DRAP representation!"]];
Return[LagkJointMomentsFromDMRAP[{H0,H1},K,L][[1]]])


LagkJointMomentsFromDMRAP[H_,K_: 0,L_:1]:=
Module[{k,M,iH0,H0p,Pl,pi,sumH,Nm,Nmm,row1,col1,mid},
	If[BuTools`CheckInput && !CheckDMRAPRepresentation[H],Throw["LagkJointMomentsFromDMRAP: input is not a valid DMRAP representation!"]];
	If[K==0,k=Length[H[[1]]]-1,k=K];
	M = Length[H]-1;
	sumH=Sum[H[[i]],{i,2,M+1}];
	iH0=Inverse[IdentityMatrix[Dimensions[H[[1]]][[1]]]-H[[1]]];
	pi=DRPSolve[iH0.sumH];
	H0p=Table[i! MatrixPower[iH0,i],{i,0,k}];
	Do[H0p[[i]]=H0p[[i]].MatrixPower[H[[1]],i-2],{i,2,Length[H0p]}];
	Pl=MatrixPower[iH0.sumH,L-1];
	Nm = {};
	Do[
		Nmm=Table[Total[pi.H0p[[i]].iH0.(H[[m+1]]).Pl.H0p[[j]]],{i,k+1},{j,k+1}];
        row1 = MomsFromFactorialMoms[Nmm[[1,2;;]]];
        col1 = MomsFromFactorialMoms[Nmm[[2;;,1]]];
        mid = JMomsFromJFactorialMoms[Nmm[[2;;,2;;]]];
        AppendTo[Nm,ArrayFlatten[{{{{Nmm[[1,1]]}}, {row1}},{Transpose[{col1}], mid}}]];
	,{m,M}];
	Return[Nm];
]


LagkJointMomentsFromDMAP[D0_,D1_,K_:0,L_:1]:=
(If[BuTools`CheckInput && !CheckDMAPRepresentation[D0,D1],Throw["LagkJointMomentsFromDMAP: Input is not a valid DMAP representation!"]];
Return[LagkJointMomentsFromDMRAP[{D0,D1},K,L][[1]]])


LagkJointMomentsFromDMMAP[D_,K_:0,L_:1]:=
(If[BuTools`CheckInput && !CheckDMMAPRepresentation[D],Throw["LagkJointMomentsFromDMMAP: input is not a valid DMMAP representation!"]];
Return[LagkJointMomentsFromDMRAP[D,K,L]])


DMAP2FromMoments[moms_,corr1_] :=
Module[ {Nm,H0,H1,oCheckInput,D0,D1},
    Nm = {{1, moms[[1]]},{moms[[1]], corr1 (moms[[2]]-moms[[1]]^2)+moms[[1]]^2}};
    {H0, H1} = DRAPFromMoments[moms, Nm];
    oCheckInput = BuTools`CheckInput;
	BuTools`CheckInput = False;
    {D0, D1} = CanonicalFromDMAP2[H0, H1];
    BuTools`CheckInput = oCheckInput;
	Return[{D0,D1}];
];


DMAPFromDRAP[H0_,H1_,prec_:N[10^-14]]:= (
If[BuTools`CheckInput && !CheckDRAPRepresentation[H0,H1],Throw["DMAPFromDRAP: Input is not a valid DRAP representation!"]];
Return[DMMAPFromDMRAP[{H0,H1},prec]];)


DMMAPFromDMRAP[H_,prec_:N[10^-14]]:= 
Module[{Transfun,Evalfun,nrep},
    Transfun[orep_, B_]:=
		Map[Inverse[B].#.B &, orep];
        
    Evalfun[orep_, k_:0]:=Module[{oH0},
        If[Mod[k,2] == 0,
            Return[-Min[orep,1-orep]];
        ,
            Return[-Total[Map[Total[Select[Flatten[#],Negative]]&,orep]]-Total[Map[Total[Select[Flatten[1-#],Negative]]&,orep]]];
        ];
	];
   
    If[BuTools`CheckInput && Not[CheckDMRAPRepresentation[H]],Throw["DMMAPFromDMRAP: Input is not a valid DMRAP representation!"]];
    nrep = FindMarkovianRepresentation[H, Transfun, Evalfun, prec];
    Return[nrep];
];


MarginalDistributionFromDRAP[H0_,H1_]:=(
If[BuTools`CheckInput && !CheckDRAPRepresentation[H0,H1],Throw["MarginalDistributionFromDRAP: Input is not a valid DRAP representation!"]];
Return[{DRPSolve[Inverse[IdentityMatrix[Dimensions[H0][[1]]]-H0].H1],H0}];)


MarginalDistributionFromDMRAP[H_]:=
Module[{hk},
If[BuTools`CheckInput && !CheckDMRAPRepresentation[H],Throw["MarginalDistributionFromDMRAP: Input is not a valid DMRAP representation!"]];
hk=Sum[H[[i]],{i,2,Length[H]}];
Return[{DRPSolve[Inverse[IdentityMatrix[Dimensions[H[[1]]][[1]]]-H[[1]]].hk],H[[1]]}];
];


MarginalDistributionFromDMAP[D0_,D1_]:=(
If[BuTools`CheckInput && !CheckDMAPRepresentation[D0,D1],Throw["MarginalDistributionFromDMAP: Input is not a valid DMAP representation!"]];
Return[MarginalDistributionFromDRAP[D0,D1]];)


MarginalDistributionFromDMMAP[D_]:=
(If[BuTools`CheckInput && !CheckDMMAPRepresentation[D],Throw["MarginalDistributionFromDMMAP: Input is not a valid DMMAP representation!"]];
Return[MarginalDistributionFromDMRAP[D]];)


MarginalMomentsFromDRAP[H0_,H1_,K_:0]:=
Module[{\[Alpha],A},
If[BuTools`CheckInput && Not[CheckDRAPRepresentation[H0,H1]], Throw["MarginalMomentsFromDRAP: Input is not a valid DRAP representation!"]];
{\[Alpha],A}=MarginalDistributionFromDRAP[H0,H1];
Return[MomentsFromMG[\[Alpha],A,K]];
]


MarginalMomentsFromDMRAP[H_,K_:0]:=
(If[BuTools`CheckInput && !CheckDMRAPRepresentation[H],Throw["MarginalMomentsFromDMRAP: Input is not a valid DMRAP representation!"]];
Return[MarginalMomentsFromDRAP[H[[1]],Sum[H[[i]],{i,2,Length[H]}],K]];)


MarginalMomentsFromDMAP[D0_,D1_,K_:0]:=
Module[{\[Alpha],d0i},
If[BuTools`CheckInput && Not[CheckDMAPRepresentation[D0,D1]], Throw["MarginalMomentsFromDMAP: Input is not a valid DMAP representation!"]];
Return[MarginalMomentsFromDRAP[D0,D1,K]];
]


MarginalMomentsFromDMMAP[D_,K_:0]:=(
If[BuTools`CheckInput && !CheckDMMAPRepresentation[D],Throw["MarginalMomentsFromDMMAP: Input is not a valid DMMAP representation!"]];
Return[MarginalMomentsFromDRAP[D[[1]],Sum[D[[i]],{i,2,Length[D]}],K]];)


DRAPFromMoments[moms_,Nm_]:=DMRAPFromMoments[moms,{Nm}];


DMRAPFromMoments[moms_,Nm_]:=
Module[ {H0,v,H0i,size,gammav,gamma1,gammavi,gamma1i,Nmi,H,H0ip},
	{v,H0}=MGFromMoments[moms];
	size=Dimensions[H0][[1]];
	H0i=Inverse[IdentityMatrix[size]-H0];
	H0ip=Prepend[Table[i! MatrixPower[H0i,i].MatrixPower[H0,i-1],{i,1,size-1}],IdentityMatrix[size]];
	gamma1=Transpose[Table[Total[H0ip[[i]],{2}],{i,size}]];
	gammav=Table[v.H0ip[[i]],{i,size}];
	gamma1i=Inverse[gamma1];
	gammavi=Inverse[gammav];
	H={H0};	
	Do[
		Nmi=Nm[[i-1]];
		Nmi=ArrayFlatten[{{{{Nmi[[1,1]]}},{FactorialMomsFromMoms[Nmi[[1,2;;]]]}},{Transpose[{FactorialMomsFromMoms[Nmi[[2;;,1]]]}], JFactorialMomsFromJMoms[Nmi[[2;;,2;;]]]}}];
		AppendTo[H,(IdentityMatrix[size]-H0).gammavi.Nmi.gamma1i];
	,{i,2,Length[Nm]+1}];
	Return[H];
];


RandomDMAP[order_, mean_: 10.,zeroEntries_: 0, maxTrials_: 1000, prec_: N[10^-7]]:=
Return[RandomDMMAP[order,1,mean,zeroEntries,maxTrials,prec]];


RandomDMMAP[order_,types_,mean_:10.,zeroEntries_:0,maxTrials_:1000, prec_: N[10^-7]]:=
Module[{trials,zeros,B,numZeros,zeroInRow,idx,a,aRowSize,d,\[Pi],fullZero,actualZeros,dd,dx,Dx,Dv,m},
If[zeroEntries>(order+1)(order-1)+types(order^2-1),Throw["RandomDMAP/DMMAP: You have given too many zeros! Try to decrease the zero entries number!"]];
actualZeros=zeroEntries;
While[actualZeros>=0,
trials=1;
While[trials<maxTrials,
zeros=Table[0,{i,order}];
numZeros=0;
idx=1;
While[numZeros<actualZeros && idx<order,
zeroInRow=RandomInteger[Min[actualZeros-numZeros,(types+1)order-1]];
zeros[[idx]]=zeroInRow;
numZeros=numZeros+zeroInRow;
idx++;
];
zeros[[order]]=actualZeros-numZeros;
If[zeros[[order]]>=(types+1)order-1,trials++;Continue[]];

B=Table[0,{i,order},{j,(types+1)order}];
aRowSize=(types+1)order;
For[idx=1,idx <= order, idx++,
a=Table[Random[],{i,aRowSize-zeros[[idx]]}];
a=a/Total[a];
a=Join[a,Table[0,{i,zeros[[idx]]}]];
a=a[[RandomSample[Range[aRowSize],aRowSize]]];(*This permutates the elements*)
B[[idx,;;]]=a[[;;]];
];(*For*)
d=Table[B[[;;,i order+1;;(i+1)order]],{i,0,types}];
(*Check the obtained representation*)

If[MatrixRank[d[[1]]-IdentityMatrix[order]] != order,trials++;Continue[]];
If[MatrixRank[Sum[d[[i]],{i,types+1}]-IdentityMatrix[order]]!= order-1,trials++;Continue[]];
\[Pi]=CTMCSolve[Sum[d[[i]],{i,types+1}]-IdentityMatrix[order]];
If[Min[Abs[\[Pi]]]<\[Epsilon],trials++;Continue[]];
fullZero=False;
Do[If[And@@Flatten[Table[Abs[d[[i,k,l]]]<\[Epsilon],{k,order},{l,order}]],fullZero=True],{i,types+1}];
If[fullZero,trials++;Continue[]];

dd = Table[Random[Real],{order}];
(* scale to the mean value *)
Dv = Map[DiagonalMatrix[1-dd].#&,d];
Dv[[1]] += DiagonalMatrix[dd];
m = MarginalMomentsFromDMMAP[Dv,1];
dx = 1 - (1-dd) m[[1]]/mean;
Dx = Map[DiagonalMatrix[1-dx].#&,d];
Dx[[1]] += DiagonalMatrix[dx];
If[Not[CheckDMMAPRepresentation[Dx]],Continue[]];
If[zeroEntries>actualZeros,
Print["RandomDMAP/DMMAP: Number of zero entries are different! Given:",zeroEntries,", in the returned representation:",actualZeros]];
Return[Dx];
];(*While*)
--actualZeros;
];(*While*)
];(*Module*)


SamplesFromDMAP[D0_,D1_,k_,initial_:Null]:=(
	If[BuTools`CheckInput && !CheckDMAPRepresentation[D0,D1],Throw["SamplesFromDMAP: Input is not a valid DMAP representation!"]];
	Return[SamplesFromDMMAP[{D0,D1},k,initial]];)


SamplesFromDMMAP[D_,k_,initial_:Null]:=
Module[{NN,cummInitial,sojourn, nextpr,r,state,stst,logp,genSamples},
	If[BuTools`CheckInput && Not[CheckDMMAPRepresentation[D]],Throw["SamplesFromDMMAP: input is not a valid DMMAP representation!"]];
    
    NN = Dimensions[D[[1]]][[1]];
	
	If[initial===Null,
		stst=MarginalDistributionFromDMMAP[D][[1]];
        cummInitial = Accumulate[stst];
        r = RandomReal[];
        state = 1;
        While[cummInitial[[state]]<=r, state++];
    ,
		state = initial;
	];
	(* auxilary variables*)
    sojourn = 1/(1-Diagonal[D[[1]]]);
	logp = Log[Diagonal[D[[1]]]];
    nextpr = DiagonalMatrix[sojourn].D[[1]];
    nextpr = nextpr - DiagonalMatrix[Diagonal[nextpr]];
    Do[nextpr = Join[nextpr, DiagonalMatrix[sojourn].D[[i]], 2],{i,2,Length[D]}];
    nextpr = Transpose[Accumulate[Transpose[nextpr]]];
    
    genSamples=Compile[{{NN,_Integer},{logp,_Real,1},{sojourn,_Real,1},{nextpr,_Real,2}},
	Module[{time,rr,cstate,nstate,x},
		x = Table[0.,{k},{2}];
		cstate = state;
		Do[
			time = 0.;
			(* play state transitions *)
			While[cstate<=NN,
				time += 1 + Floor[Log[RandomReal[]] / logp[[cstate]]];
				rr = RandomReal[];
				nstate = 1;
				While[nextpr[[cstate,nstate]]<=rr, nstate++];
				cstate = nstate;
			];
			x[[n,1]]=time; 
			x[[n,2]]=Ceiling[cstate/NN]-1;
			cstate=Mod[cstate-1,NN]+1;
		,{n,k}];
		Return[x]
	]];
	If[Length[D]>2,
		Return[genSamples[NN,logp,sojourn,nextpr]]
	, 
		Return[genSamples[NN,logp,sojourn,nextpr][[;;,1]]]
	];
];


End[(* Private *)];
EndPackage[];
