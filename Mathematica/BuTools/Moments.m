(* ::Package:: *)

(*
   BuTools Moments Package
*)

BeginPackage["BuTools`Moments`"];
NormMomsFromMoms::usage = "nm = NormMomsFromMoms[m]: Returns the normalized moments given the raw moments";
MomsFromNormMoms::usage = "m = MomsFromNormMoms[nm]: Returns the raw moments given the normalized moments";
ReducedMomsFromMoms::usage = "rm = ReducedMomsFromMoms[m]: Returns the reduced moments given the raw moments";
MomsFromReducedMoms::usage = "m = MomsFromReducedMoms[rm]: Returns the raw moments given the reduced moments";
FactorialMomsFromMoms::usage = "fm = FactorialMomsFromMoms[m]: Returns the factorial moments given the raw moments";
MomsFromFactorialMoms::usage = "m = MomsFromFactorialMoms[fm]: Returns the raw moments given the factorial moments";
HankelMomsFromMoms::usage = "hm = HankelMomsFromMoms[m]: Returns the Hankel moments given the raw moments";
MomsFromHankelMoms::usage = "m = MomsFromHankelMoms[hm]: Returns the raw moments given the Hankel moments";
JFactorialMomsFromJMoms::usage = "jfm = JFactorialMomsFromJMoms[jm]: Returns the joint factorial moments given the joint raw moments";
JMomsFromJFactorialMoms::usage = "jm = JMomsFromJFactorialMoms[jfm]: Returns the joint raw moments given the joint factorial moments";
CheckMoments::usage = "r = CheckMoments[m, prec]: Checks if the given moment sequence belongs to a distribution with support (0,inf)";


Begin["`Private`"];


If[Not[ValueQ[BuTools`CheckInput]],BuTools`CheckInput=True];
If[Not[ValueQ[BuTools`CheckPrecision]],BuTools`CheckPrecision=N[10^-12]];
If[Not[ValueQ[BuTools`Verbose]],BuTools`Verbose=False];


NormMomsFromMoms[ moms_]:=Prepend[Table[moms[[i]]/(moms[[i-1]] moms[[1]]),{i,2,Length[moms]}],moms[[1]]];


MomsFromNormMoms[ normmoms_] :=
Module[ {moms},
moms={normmoms[[1]]};
Do[AppendTo[moms,normmoms[[i]] moms[[i-1]] moms[[1]]],{i,2,Length[normmoms]}];
Return[moms];
];


ReducedMomsFromMoms[moms_]:= Table[moms[[i]]/i!,{i,Length[moms]}];


MomsFromReducedMoms[ rmoms_] := Table[rmoms[[i]] i!,{i,Length[rmoms]}];


FactorialMomsFromMoms[moms_]:=
Table[Table[Coefficient[Product[(x-k),{k,0,i-1}],x,j],{j,i}].moms[[1;;i]],{i,Length[moms]}];
(*Module[{fmoms,eh,i},
fmoms={};
Do[
eh=Table[Coefficient[Product[(x-k),{k,0,i-1}],x,j],{j,i}];
AppendTo[fmoms,eh.moms[[1;;i]]];
,{i,Length[moms]}];
Return[fmoms];
];*)


MomsFromFactorialMoms[fmoms_]:=
Module[{moms,eh},
moms={fmoms[[1]]};
Do[
AppendTo[moms,fmoms[[i]]-Table[Coefficient[Product[(x-k),{k,0,i-1}],x,j],{j,i-1}].moms[[1;;i-1]]];
,{i,2,Length[fmoms]}];
Return[moms];
];


HankelMomsFromMoms[moms_]:=
Module[{hm,NN,H,i},
hm={};
Do[
	If[Mod[i,2]==1,
		NN=(i-1)/2+1;
		H=HankelMatrix[moms[[1;;NN]],moms[[NN;;2 NN-1]]];
	,
		NN=i/2+1;
		H=HankelMatrix[Prepend[moms[[1;;NN-1]],1],moms[[NN-1;;2 NN-2]]];
	];
	AppendTo[hm,Det[H]];
,{i,Length[moms]}];
Return[hm];
];


MomsFromHankelMoms[hmoms_]:=
Module[{moms,NN,H,i,j,h,rH,rHd,cofactor},
moms={hmoms[[1]]};
Do[
	If[Mod[i,2]==0,
		NN=i/2+1;
		H=HankelMatrix[moms[[1;;NN]],Append[moms[[NN;;2 NN-2]],0]];
	,
		NN=(i+1)/2+1;
		H=HankelMatrix[Prepend[moms[[1;;NN-1]],1],Append[moms[[NN-1;;2 NN-3]],0]];
	];
	h=hmoms[[i+1]];
	rH = H[[1;;NN-1,;;]];
	Do[
		rHd = rH;
		rHd=Drop[rHd,{},{j+1}];
		cofactor = (-1)^(NN+j-1) Det[rHd];
		If[j<NN-1,h=h-cofactor H[[NN,j+1]],AppendTo[moms,h/cofactor]];
	,{j,0,NN-1}];
,{i,Length[hmoms]-1}];
Return[moms];
];


JFactorialMomsFromJMoms[jmoms_]:=
Module[{size1,size2,i,j,prod,eh,ret},
size1=Dimensions[jmoms][[1]];
size2=Dimensions[jmoms][[2]];
ret=Table[0,{k,1,size1},{l,1,size2}];
For[i=1,i<=size1,++i,
	For[j=1,j<=size2,++j,
		prod=Product[(x-k),{k,0,i-1}]Product[(y-k),{k,0,j-1}];
		eh=Table[Coefficient[prod,x^k y^l],{k,1,i},{l,1,j}];
		ret[[i,j]]=Tr[jmoms[[1;;i,1;;j]].Transpose[eh]];
	];
];
Return[ret];
];


JMomsFromJFactorialMoms[jfmoms_]:=
Module[{size1,size2,i,j,prod,eh,ret},
size1=Dimensions[jfmoms][[1]];
size2=Dimensions[jfmoms][[2]];
ret=Table[0,{k,size1},{l,size2}];
For[i=1,i<=size1,++i,
	For[j=1,j<=size2,++j,
		prod=Product[(x-k),{k,0,i-1}]Product[(y-k),{k,0,j-1}];
		eh=-Table[Coefficient[prod,x^k y^l],{k,1,i},{l,1,j}];
		ret[[i,j]]=jfmoms[[i,j]]+Tr[ret[[1;;i,1;;j]].Transpose[eh]];
	];
];
Return[ret];
];


CheckMoments[moms_,precision_:Null]:=
Module[{NN,H,H0,m,res,prec},
If[Not[NumericQ[precision]], prec=BuTools`CheckPrecision, prec=precision];
If[BuTools`CheckInput && Mod[Length[moms],2]==0,Throw["CheckMoments: the number of moments must be odd!"]];
m=Prepend[moms,1];
NN=Floor[Length[m]/2]-1;
res = True;
Do[
	H=HankelMatrix[m[[1;;n+1]],m[[n+1;;2 n+1]]];
	H0=HankelMatrix[m[[2;;n+2]],m[[n+2;;2 n+2]]];
	If[Det[H]<-prec || Det[H0]<-prec, res=False];
,{n,0,NN}];
Return[res];
];


End[(* Private *)];

EndPackage[];
