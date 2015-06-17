(* ::Package:: *)

(*
   BuTools_MAP Package
*)


BeginPackage["BuTools`MAP`"];
CanonicalFromMAP2::usage = "alma";
CheckMAPRepresentation::usage = "alma";
CheckMMAPRepresentation::usage = "alma";
CheckRAPRepresentation::usage = "alma";
CheckMRAPRepresentation::usage = "alma";
LagCorrelationsFromRAP::usage = "alma";
LagCorrelationsFromMAP::usage = "alma";
LagkJointMomentsFromRAP::usage = "alma";
LagkJointMomentsFromMRAP::usage = "alma";
LagkJointMomentsFromMAP::usage = "alma";
LagkJointMomentsFromMMAP::usage = "alma";
MAP2CorrelationBounds::usage = "alma";
MAP2FromMoments::usage = "alma";
MAPFromFewMomentsAndCorrelations::usage = "alma";
MAPFromRAP::usage = "alma";
MMAPFromMRAP::usage = "alma";
MarginalDistributionFromRAP::usage = "alma";
MarginalDistributionFromMRAP::usage = "alma";
MarginalDistributionFromMAP::usage = "alma";
MarginalDistributionFromMMAP::usage = "alma";
MarginalMomentsFromRAP::usage = "alma";
MarginalMomentsFromMRAP::usage = "alma";
MarginalMomentsFromMAP::usage = "alma";
MarginalMomentsFromMMAP::usage = "alma";
MinimalRepFromRAP::usage = "alma";
MinimalRepFromMRAP::usage = "alma";
RAPFromMoments::usage = "alma";
MRAPFromMoments::usage = "alma";
RandomMAP::usage = "alma";
RandomMMAP::usage = "alma";
RAPFromMomentsAndCorrelations::usage = "alma";
SamplesFromMAP::usage = "alma";
SamplesFromMMAP::usage = "alma";


Begin["`Private`"];


Needs["BuTools`RepTrans`"];
Needs["BuTools`Moments`"];
Needs["BuTools`MC`"];
Needs["BuTools`PH`"];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckInput"]],BuTools`CheckInput=True];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckPrecision"]],BuTools`CheckPrecision=N[10^-12]];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`Verbose"]],BuTools`Verbose=False];


CanonicalFromMAP2[D0_,D1_]:=
Module[{moms,corr1},
	If[BuTools`CheckInput && Dimensions[D0][[1]]!=2, Throw["CanonicalFromMAP2: Size is not 2!"]];
	If[BuTools`CheckInput && Not[CheckMAPRepresentation[D0, D1]],Throw["CanonicalFromMAP2: Input is not a valid MAP representation!"]];
    moms = MarginalMomentsFromMAP[D0, D1, 3];
    corr1 = LagCorrelationsFromMAP[D0, D1, 1];
    Return[MAP2FromMoments[moms, corr1]];
];


CheckMAPRepresentation[D0_,D1_,precision_:Null] :=
Module[{prec},
	If[Not[NumericQ[precision]], prec=BuTools`CheckPrecision, prec=precision];
	If[Not[CheckGenerator[D0,True,prec]],Return[False]];
	If[Dimensions[D0]!=Dimensions[D1],
		If[BuTools`Verbose,Print["CheckMAPRepresentation: D0 and D1 have different sizes!"]]; 
		Return[False]
	];
	If[ Min[D1]<-prec,
		If[BuTools`Verbose,Print["CheckMAPRepresentation: D1 has negative element!"]];
		Return[False]
	];
	If[Max[Abs[Total[D0+D1,{2}]]]>prec,
		If[BuTools`Verbose==True,Print["CheckMAPRepresentation: A rowsum of D0+D1 is not 0 (at precision ",prec,")!"]];
		Return[False]
	];
	Return[True]
];


CheckMMAPRepresentation[D_,precision_:Null]:=
Module[{prec},
	If[Not[NumericQ[precision]], prec=BuTools`CheckPrecision, prec=precision];
	If[Min[D[[2;;]]]<-prec,
		If[BuTools`Verbose,Print["CheckMMAPRepresentation: One of the matrices has a negative element!"]];
		Return[False]
	];
	Return[CheckMAPRepresentation[D[[1]],Sum[D[[i]],{i,2,Length[D]}],prec]];
];


CheckRAPRepresentation[D0_,D1_, precision_:Null]:=
Module[ {h,size1,size2,size,eig,ceig,reig,prec},
	If[Not[NumericQ[precision]], prec=BuTools`CheckPrecision, prec=precision];
	If[ Dimensions[D0][[1]]!=Dimensions[D0][[2]],
		If[BuTools`Verbose,Print["CheckRAPRepresentation: D0 is not a quadratic matrix!"]];
		Return[False]
	];
	If[Dimensions[D1][[1]]!=Dimensions[D1][[2]],
		If[BuTools`Verbose,Print["CheckRAPRepresentation: D1 is not a quadratic matrix!"]];
		Return[False]
	];
	If[Dimensions[D0]!=Dimensions[D1],
		If[BuTools`Verbose,Print["CheckRAPRepresentation: D0 and D1 have different sizes!"]]; 
		Return[False]
	];
	If[Max[Abs[Total[D0+D1,{2}]]] > prec,
		If[BuTools`Verbose,Print["CheckRAPRepresentation: A rowsum of D0+D1 is not 0! (precision:",prec,")"]];
		Return[False]
	];
	eig=Eigenvalues[D0//N];
	If[Max[Re[eig]]>=-prec,
	If[BuTools`Verbose,Print["CheckRAPRepresentation: there is an eigenvalue of D0 with non-negative real part (at precision ",prec,")!"]];
		Return[False]
	];
	size=Dimensions[D0][[1]];
	ceig=Table[Switch[Element[eig[[i]],Reals],True,-Infinity,False,Re[eig[[i]]]],{i,size}];
	reig=Table[Switch[Element[eig[[i]],Reals],True,Re[eig[[i]]],False,-Infinity],{i,size}];

	If[Max[reig]<Max[ceig],
		If[BuTools`Verbose,Print["CheckRAPRepresentation: The dominant eigenvalue of D0 is not real!"]];
		Return[False]
	];
	If[Max[reig]==Max[ceig] && BuTools`Verbose,
		Print["CheckRAPRepresentation: The dominant and a complex eigenvalue of matrix0 has the same real part!"]
	];
	Return[True]
];


CheckMRAPRepresentation[H_, precision_:Null] := 
	CheckRAPRepresentation[H[[1]],Sum[H[[i]],{i,2,Length[H]}],precision];


LagCorrelationsFromRAP[H0_,H1_,L_:1]:=
Module[{H0i,P,pi,moms,acf},
	If[BuTools`CheckInput && !CheckRAPRepresentation[H0,H1],Throw["LagCorrelationsFromRAP: Input is not a valid RAP representation!"]];
	H0i = Inverse[-H0];
	P = H0i.H1;
	pi = DRPSolve[P];
	moms = MomentsFromME[pi,H0];
	pi = pi.H0i.P.
	acf = {};
	Do[
		AppendTo[acf,(Total[pi.H0i]-moms[[1]]^2) / (moms[[2]]-moms[[1]]^2)];
		pi = pi.P;
	,{L}];
	Return[acf];
];


LagCorrelationsFromMAP[D0_,D1_,L_:1]:=
(If[BuTools`CheckInput && !CheckMAPRepresentation[D0,D1],Throw["LagCorrelationsFromMAP: Input is not a valid MAP representation!"]];
Return[LagCorrelationsFromRAP[D0,D1,L]])


LagkJointMomentsFromRAP[H0_,H1_,K_:0, L_:1]:=
(If[BuTools`CheckInput && !CheckRAPRepresentation[H0,H1],Throw["LagkJointMomentsFromRAP: Input is not a valid RAP representation!"]];
Return[LagkJointMomentsFromMRAP[{H0,H1},K,L][[1]]])


LagkJointMomentsFromMRAP[H_,K_: 0,L_:1]:=
Module[{k,M,iH0,H0p,Pl,pi,sumH,Nm},
	If[BuTools`CheckInput && !CheckMRAPRepresentation[H],Throw["LagkJointMomentsFromMRAP: input is not a valid MRAP representation!"]];
	If[K==0,k=Length[H[[1]]]-1,k=K];
	M = Length[H]-1;
	sumH=Sum[H[[i]],{i,2,M+1}];
	iH0=Inverse[-H[[1]]];
	pi=DRPSolve[iH0.sumH];
	H0p=Table[i! MatrixPower[iH0,i],{i,0,k}];
	Pl=MatrixPower[iH0.sumH,L-1];
	Nm=Table[Table[Total[pi.H0p[[i]].iH0.(H[[m+1]]).Pl.H0p[[j]]],{i,k+1},{j,k+1}],{m,M}];
	Return[Nm];
]


LagkJointMomentsFromMAP[D0_,D1_,K_:0,L_:1]:=
(If[BuTools`CheckInput && !CheckMAPRepresentation[D0,D1],Throw["LagkJointMomentsFromMAP: Input is not a valid MAP representation!"]];
Return[LagkJointMomentsFromMRAP[{D0,D1},K,L][[1]]])


LagkJointMomentsFromMMAP[D_,K_:0,L_:1]:=
(If[BuTools`CheckInput && !CheckMMAPRepresentation[D],Throw["LagkJointMomentsFromMMAP: input is not a valid MMAP representation!"]];
Return[LagkJointMomentsFromMRAP[D,K,L]])


MAP2CorrelationBounds[moms_]:=
Module[{m1,m2,m3,h2,h3,cv2,gub,glb},
    {m1, m2, m3} = moms;
    h2 = m2/(2 m1 m1) - 1;
    h3 = m3/(6 m1 m1 m1)-m2 m2/(4 m1 m1 m1 m1);
    cv2 = m2/m1/m1 - 1;
    If[h2>=0, gub=h2, gub=-(h2+Sqrt[-h3])^2];
    If[h2<=0 || h3/h2+h2<1, glb=-h3-h2 h2, glb=h2 (h3+h2 h2-h2-Sqrt[(h3+h2 h2-h2)^2+4 h2 h2 h2]) / (h3+h2 h2-h2+Sqrt[(h3+h2 h2-h2)^2+4 h2 h2 h2])];
    If[h2>=0, Return[{glb/cv2,gub/cv2}], Return[{gub/cv2,glb/cv2}]];
];


MAP2FromMoments[moms_,corr1_] :=
Module[ {m1,m2,m3,\[Lambda]1,\[Lambda]2,p,\[Gamma],tau,T,\[Alpha],corrl,corru,a,b},
	If[Abs[m2-2 m1 m1] < BuTools`CheckPrecision && Abs[corr1] > BuTools`CheckPrecision, Throw["We do not allow correlation in case of exponentially distributed marginal"]];

	{m1,m2,m3}=moms;
	(* Perform PH fitting *)
	{tau,T} = PH2From3Moments[moms];
	\[Lambda]1 = -T[[1,1]];
	\[Lambda]2 = -T[[2,2]];
	p = tau[[1]];
	\[Alpha] = \[Lambda]1/\[Lambda]2;
    (* Check the feasibility of the correlation parameter*)
    {corrl, corru} = MAP2CorrelationBounds[moms];
    If[corr1<corrl,Throw["The correlation parameter is too small!"]];
    If[corr1>corru,Throw["The correlation parameter is too large!"]];
    \[Gamma] = corr1 (m2-m1 m1) / (m2/2 - m1 m1);
	(* Perform matching *)
	If[\[Gamma]>= 0,
		a=1/(2 \[Alpha]) (1+\[Alpha] \[Gamma] -p (1-\[Gamma])-Sqrt[(1+\[Alpha] \[Gamma] -p (1-\[Gamma]))^2-4 \[Alpha] \[Gamma]]);
		b=1/2  (1+\[Alpha] \[Gamma] -p (1-\[Gamma])+Sqrt[(1+\[Alpha] \[Gamma] -p (1-\[Gamma]))^2-4 \[Alpha] \[Gamma]]);
		Return[{{{-\[Lambda]1, (1-a)\[Lambda]1},{0, -\[Lambda]2}}, {{a \[Lambda]1, 0},{(1-b) \[Lambda]2, b \[Lambda]2}}}];
	];
	If[ \[Gamma]<0,
		a=-\[Gamma]/(p (1-\[Gamma])- \[Alpha] \[Gamma]);
		b=p (1-\[Gamma])- \[Alpha] \[Gamma];
		Return[{{{-\[Lambda]1, (1-a)\[Lambda]1},{0, -\[Lambda]2}}, {{0, a \[Lambda]1},{b \[Lambda]2, (1-b) \[Lambda]2}}}];
	];
];


MAPFromFewMomentsAndCorrelations[moms_,corr1_,rr_:None]:=
Module[{m1,c2,l3,\[Rho]1,r,m11,m12,p1,p2,cv21,cv22,l31,l32,\[Alpha]1,A1,\[Alpha]2,A2,N1,N2,h1,h2,D0,D1,m21,m22,m31,m32},
m1 = moms[[1]];
c2 = moms[[2]]/moms[[1]]^2 - 1;
If[Length[moms]>2, l3=moms[[3]] moms[[1]]/moms[[2]]^2 - 1, l3 = None];   
\[Rho]1=corr1;
If[\[Rho]1>=0,
If[rr=!=None && rr<=0 && rr>=1,Print["Parameter r is out of range!"];Abort[]];
If[rr===None, 
r=(2 \[Rho]1)/(1+\[Rho]1);
p1= 1/(1+c2) (1-(1+\[Rho]1)/2);
p2= c2/(1+c2) (1-(1+\[Rho]1)/2);
,
r=rr;
p1= 1/(1+c2) (1-\[Rho]1/r);
p2= c2/(1+c2) (1-\[Rho]1/r);
];
m11 = m1(1-Sqrt[r]);
m12 = m1(1+c2 Sqrt[r]);
If[l3===None,
cv21=(Sqrt[c2] (1+c2) (1+Sqrt[r]))/(1-Sqrt[c2] (-1+Sqrt[r])+c2 Sqrt[r]);
cv22=-((c2 (1+c2) (-1+r))/((1+c2 Sqrt[r]) (1-Sqrt[c2] (-1+Sqrt[r])+c2 Sqrt[r])));
(*Print[r,", ",m11,", ", m12, ", ", cv21, ", ",cv22];*)
m21 = (cv21+1) m11^2;
m22 = (cv22+1) m12^2;
{\[Alpha]1,A1}=APHFrom2Moments[{m11,m21}];
{\[Alpha]2,A2}=APHFrom2Moments[{m12,m22}];
,
cv21=(c2+Sqrt[r])/(1-Sqrt[r]);
cv22=c2 (1-Sqrt[r])/(1+c2 Sqrt[r]);
l31=((1+c2) l3)/(c2(1-Sqrt[r])+Sqrt[(1 +c2 Sqrt[r])c2(1-Sqrt[r])]);
l32=((1+c2) l3)/((1+c2 Sqrt[r])+Sqrt[(1 +c2 Sqrt[r])c2(1-Sqrt[r])]);
(*Print[r,", ",m11,", ", m12, ", ", cv21, ", ",cv22,", ",l31,", ", l32];*)
If[cv21<=0||cv22<=0||l31<=0 || l32<=0, Print["Not feasible. c2=",c2,", l3=",l3, ", r=",r]; Return[{{},{}}]];
m21 = (cv21+1)m11^2;
m22 = (cv22+1)m12^2;
m31 = (l31+1)m21^2/m11;
m32 = (l32+1)m22^2/m12;
{\[Alpha]1,A1}=APHFrom3Moments[{m11,m21,m31}];
{\[Alpha]2,A2}=APHFrom3Moments[{m12,m22,m32}];
],
If[c2>=1,
If[rr=!=None && rr<=0 && rr>=1/c2,Print["Parameter r is out of range!"];Abort[]];
If[rr===None, 
r=-((2 \[Rho]1)/(1-c2 \[Rho]1));
p1= 1/2 (1+(1-c2 \[Rho]1)/2);
p2= 1/2 (1+(1-c2 \[Rho]1)/2);
,
r=rr;
p1= 1/2 (1-\[Rho]1/r);
p2= 1/2 (1-\[Rho]1/r);
];
,
If[rr=!=None && rr<=0 && rr>=1,Print["Parameter r is out of range!"];Abort[]];
If[rr===None, 
r=-((2 \[Rho]1)/(1-\[Rho]1));
p1= 1/2 (1+(1-\[Rho]1)/2);
p2= 1/2 (1+(1-\[Rho]1)/2);
,
r=rr;
p1= 1/2 (1-\[Rho]1/r);
p2= 1/2 (1-\[Rho]1/r);];
];
m11 = m1(1-Sqrt[c2 r]);
m12 = m1(1+ Sqrt[c2 r]);
If[l3===None,
cv21=c2 (1-r)/(1-Sqrt[c2 r]);
cv22=c2 (1-r)/(1+Sqrt[c2 r]);
(*Print[r,", ",m11,", ", m12, ", ", cv21, ", ",cv22];*)
m21 = (cv21+1) m11^2;
m22 = (cv22+1) m12^2;
{\[Alpha]1,A1}=APHFrom2Moments[{m11,m21}];
{\[Alpha]2,A2}=APHFrom2Moments[{m12,m22}];
,
cv21=(c2+Sqrt[c2 r])/(1-Sqrt[c2 r]);
cv22=(c2-Sqrt[c2 r])/(1+Sqrt[c2 r]);
l31=2 l3 1/(1-Sqrt[c2 r]+Sqrt[1-c2 r]);
l32=2 l3 1/(1+Sqrt[c2 r]+Sqrt[1-c2 r]);
If[cv21<=0||cv22<=0||l31<=0 || l32<=0, Print["Not feasible. c2=",c2,", l3=",l3, ", r=",r]; Return[{{},{}}]];
(*Print[r,", ",m11,", ", m12, ", ", cv21, ", ",cv22,", ",NumberForm[l31,50],", ", l32];*)
m21 = (cv21+1)m11^2;
m22 = (cv22+1)m12^2;
m31 = (l31+1)m21^2/m11;
m32 = (l32+1)m22^2/m12;
{\[Alpha]1,A1}=APHFrom3Moments[{m11,m21,m31}];
{\[Alpha]2,A2}=APHFrom3Moments[{m12,m22,m32}];
];
];
N1 = Length[\[Alpha]1];
N2 = Length[\[Alpha]2];
h1=Table[{1},{N1}];
h2=Table[{1},{N2}];
D0=Table[0,{N1+N2},{N1+N2}];
D0[[1;;N1,1;;N1]]=A1;
D0[[N1+1;;N1+N2,N1+1;;N1+N2]]=A2;
D1=Table[0,{N1+N2},{N1+N2}];
D1[[1;;N1,1;;N1]]=-A1.h1.{\[Alpha]1}(1-p1);
D1[[1;;N1,N1+1;;N1+N2]]=-A1.h1.{\[Alpha]2} p1;
D1[[N1+1;;N1+N2,1;;N1]]=-A2.h2.{\[Alpha]1} p2;
D1[[N1+1;;N1+N2,N1+1;;N1+N2]]=-A2.h2.{\[Alpha]2}(1-p2);
Return[{D0,D1}];
];


MAPFromRAP[H0_,H1_,prec_:N[10^-14]]:= (
If[BuTools`CheckInput && !CheckRAPRepresentation[H0,H1],Throw["MAPFromRAP: Input is not a valid RAP representation!"]];
Return[MMAPFromMRAP[{H0,H1},prec]];)


MMAPFromMRAP[H_,prec_:N[10^-14]]:= 
Module[{Transfun,Evalfun,nrep},
    Transfun[orep_, B_]:=
		Map[Inverse[B].#.B &, orep];
        
    Evalfun[orep_, k_:0]:=Module[{oH0},
        oH0 = orep[[1]] - DiagonalMatrix[Diagonal[orep[[1]]]];
        If[Mod[k,2] == 0,
            Return[-Min[oH0,orep[[2;;]]]];
        ,
            Return[-Total[Map[Total[Select[#,Negative]]&,orep[[2;;]]]]-Total[Select[oH0,Negative]]];
        ];
	];
   
    If[BuTools`CheckInput && Not[CheckMRAPRepresentation[H]],Throw["MMAPFromMRAP: Input is not a valid MRAP representation!"]];
    nrep = FindMarkovianRepresentation[H, Transfun, Evalfun, prec];
    Return[nrep];
];


MarginalDistributionFromRAP[H0_,H1_]:=(
If[BuTools`CheckInput && !CheckRAPRepresentation[H0,H1],Throw["MarginalDistributionFromRAP: Input is not a valid RAP representation!"]];
Return[{DRPSolve[Inverse[-H0].H1],H0}];)


MarginalDistributionFromMRAP[H_]:=
Module[{hk},
If[BuTools`CheckInput && !CheckMRAPRepresentation[H],Throw["MarginalDistributionFromMRAP: Input is not a valid MRAP representation!"]];
hk=Sum[H[[i]],{i,2,Length[H]}];
Return[{DRPSolve[Inverse[-H[[1]]].hk],H[[1]]}];
];


MarginalDistributionFromMAP[D0_,D1_]:=(
If[BuTools`CheckInput && !CheckMAPRepresentation[D0,D1],Throw["MarginalDistributionFromMAP: Input is not a valid MAP representation!"]];
Return[MarginalDistributionFromRAP[D0,D1]];)


MarginalDistributionFromMMAP[D_]:=
(If[BuTools`CheckInput && !CheckMMAPRepresentation[D],Throw["MarginalDistributionFromMMAP: Input is not a valid MMAP representation!"]];
Return[MarginalDistributionFromMRAP[D]];)


MarginalMomentsFromRAP[H0_,H1_,K_:0]:=
Module[{\[Alpha],A},
If[BuTools`CheckInput && Not[CheckRAPRepresentation[H0,H1]], Throw["MarginalMomentsFromRAP: Input is not a valid RAP representation!"]];
{\[Alpha],A}=MarginalDistributionFromRAP[H0,H1];
Return[MomentsFromME[\[Alpha],A,K]];
]


MarginalMomentsFromMRAP[H_,K_:0]:=
(If[BuTools`CheckInput && !CheckMRAPRepresentation[H],Throw["MarginalMomentsFromMRAP: Input is not a valid MRAP representation!"]];
Return[MarginalMomentsFromRAP[H[[1]],Sum[H[[i]],{i,2,Length[H]}],K]];)


MarginalMomentsFromMAP[D0_,D1_,K_:0]:=
Module[{\[Alpha],d0i},
If[BuTools`CheckInput && Not[CheckMAPRepresentation[D0,D1]], Throw["MarginalMomentsFromMAP: Input is not a valid MAP representation!"]];
Return[MarginalMomentsFromRAP[D0,D1,K]];
]


MarginalMomentsFromMMAP[D_,K_:0]:=(
If[BuTools`CheckInput && !CheckMMAPRepresentation[D],Throw["MarginalMomentsFromMMAP: Input is not a valid MMAP representation!"]];
Return[MarginalMomentsFromRAP[D[[1]],Sum[D[[i]],{i,2,Length[D]}],K]];)


MinimalRepFromMRAP[H_,how_:"obscont",precision:N[10^-12]]:= 
Module[{B,n,alpha,A,T,R},
If[BuTools`CheckInput && !CheckMRAPRepresentation[H],Throw["MinimalRepFromMRAP: Input is not a valid MRAP representation!"]];
Which[
	how=="cont",
        {B, n} = MStaircase[H, Table[1,{Dimensions[H[[1]]][[1]]},{1}], precision];
		R=Map[(Inverse[B].#.B)[[1;;n,1;;n]]&,H];,
	how=="obs",
		{alpha, A} = MarginalDistributionFromMRAP[H];
        {B, n} = MStaircase[Map[Transpose,H], Transpose[{alpha}], precision];
        R=Map[(Inverse[B].#.B)[[1;;n,1;;n]]&,H];,
	how=="obscont",
        T = MinimalRepFromMRAP[H, "cont", precision];
        R = MinimalRepFromMRAP[T, "obs", precision];
];
Return[R];
];


MinimalRepFromRAP[H0_,H1_,how_:"obscont",precision:N[10^-12]]:= 
MinimalRepFromMRAP[{H0,H1},how,precision];


RAPFromMoments[moms_,Nm_]:=MRAPFromMoments[moms,{Nm}];


MRAPFromMoments[moms_,Nm_]:=
Module[ {H0,v,H0i,size,gammav,gamma1,gammavi,gamma1i},

	{v,H0}=MEFromMoments[moms];
	H0i=Inverse[-H0];
	size=Dimensions[H0][[1]];
	gamma1=Transpose[Table[(i-1)! Total[MatrixPower[H0i,i-1],{2}],{i,size}]];
	gammav=Table[(i-1)! v.MatrixPower[H0i,i-1],{i,size}];
	gamma1i=Inverse[gamma1];
	gammavi=Inverse[gammav];	
	Return[Prepend[Map[-H0.gammavi.#.gamma1i &, Nm],H0]];
];


RandomMAP[order_, mean_: 1.,zeroEntries_: 0, maxTrials_: 1000, prec_: N[10^-7]]:=
Return[RandomMMAP[order,1,mean,zeroEntries,maxTrials,prec]];


RandomMMAP[order_,types_,mean_:1.,zeroEntries_:0,maxTrials_:1000, prec_: N[10^-7]]:=
Module[{trials,zeros,B,numZeros,zeroInRow,idx,a,rp,aRowSize,d,\[Alpha],A,mom,\[Pi],fullZero,actualZeros},
If[types<1,Throw["RandomMMAP: 'types' must be positive integer!"]];
If[zeroEntries>(order+1)(order-1)+types(order^2-1),Throw["RandomMAP/MMAP: You have given too many zeros! Try to decrease the zero entries number!"]];
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
aRowSize=(types+1)order-1;
For[idx=1,idx <= order, idx++,
a=Table[Random[],{i,aRowSize-zeros[[idx]]}];
a=Join[a,Table[0,{i,zeros[[idx]]}]];
a=a[[RandomSample[Range[aRowSize],aRowSize]]];(*This permutates the elements*)
B[[idx,1;;idx-1]]=a[[1;;idx-1]];
B[[idx,idx+1;;]]=a[[idx;;]];
];(*For*)
d=Table[B[[;;,i order+1;;(i+1)order]],{i,0,types}];
Do[d[[1,i,i]]=-Total[B[[i]]],{i,1,order}];

(*Check the obtained representation*)

If[MatrixRank[d[[1]]] != order,trials++;Continue[]];
If[MatrixRank[Sum[d[[i]],{i,types+1}]]!= order-1,trials++;Continue[]];
\[Pi]=CTMCSolve[Sum[d[[i]],{i,types+1}]];
If[Min[Abs[\[Pi]]]<prec,trials++;Continue[]];
fullZero=False;
Do[If[And@@Flatten[Table[Abs[d[[i,k,l]]]<prec,{k,order},{l,order}]],fullZero=True],{i,types+1}];
If[fullZero,trials++;Continue[]];

{\[Alpha],A}=MarginalDistributionFromMMAP[d];
mom=MomentsFromPH[\[Alpha],A][[1]];
d=d mom/mean;
If[zeroEntries>actualZeros,
Print["RandomMAP/MMAP: Number of zero entries are different! Given:",zeroEntries,", in the returned representation:",actualZeros];
];
Return[d];
];(*While*)
--actualZeros;
];(*While*)
];(*Module*)


(*calculates a RAP representation based on 2n-1 moments and 2n-3 autocorrelations
using the K.Mitchell,A.van de Liefvoort method *) 
RAPFromMomentsAndCorrelations[moms_,corr_]:=
Module[{M,D0,alpha,rcorr,tmp,X,NN, T1,T2,U1,U2,Y,II},

    {alpha, D0} = MEFromMoments[moms];
    M = Length[alpha];
	
	If[Dimensions[corr]<2 M -3, Throw["RAPFromMomentsAndCorrelations: There are less then 2n-3 correlation parameters given"]];

	rcorr=corr[[1;;2 M -3]]/((moms[[2]]/2-moms[[1]]^2)/(moms[[2]]-moms[[1]]^2));
	rcorr = MomsFromReducedMoms[rcorr];
    {tmp,X} = MEFromMoments[rcorr];
	NN = Dimensions[X][[1]];
	If[(NN+1)!= M, Throw["RAPFromMomentsAndCorrelations: Correlation order is different from ME order"]];

	T1=Table[If[j<=i,1,0],{i,NN},{j,NN}];
	T2=Table[If[j<=i,1,0],{i,M},{j,M}];
	U1=Table[If[j>=i,1/(NN-i+1),0],{i,NN},{j,NN}];
	U2=Table[If[j>=i,1/(M-i+1),0],{i,M},{j,M}];

	Y = -Inverse[T1].U1.Inverse[X].Inverse[U1].T1;
	II=IdentityMatrix[NN+1];
	II[[2;;,2;;]]=Y;
	Y = II;

	Return[{D0,-D0.Inverse[U2].T2.Y.Inverse[T2].U2}];
];


SamplesFromMAP[D0_,D1_,k_,initial_:Null]:=(
	If[BuTools`CheckInput && !CheckMAPRepresentation[D0,D1],Throw["SamplesFromMAP: Input is not a valid MAP representation!"]];
	Return[SamplesFromMMAP[{D0,D1},k,initial]];)


SamplesFromMMAP[D_,k_,initial_:Null]:=
Module[{NN,cummInitial,sojourn, nextpr,x,time,r,state,nstate},
	If[BuTools`CheckInput && Not[CheckMMAPRepresentation[D]],Throw["SamplesFromMMAP: input is not a valid MMAP representation!"]];
    
    NN = Dimensions[D[[1]]][[1]];
	
	If[initial===Null,
		stst=MarginalDistributionFromMMAP[D][[1]];
        cummInitial = Accumulate[stst];
        r = RandomReal[];
        state = 1;
        While[cummInitial[[state]]<=r, state++];
    ,
		state = initial;
	];
	(* auxilary variables*)
    sojourn = -1/Diagonal[D[[1]]];
    nextpr = DiagonalMatrix[sojourn].D[[1]];
    nextpr = nextpr - DiagonalMatrix[Diagonal[nextpr]];
    Do[nextpr = Join[nextpr, DiagonalMatrix[sojourn].D[[i]], 2],{i,2,Length[D]}];
    nextpr = Transpose[Accumulate[Transpose[nextpr]]];
    
    If[Length[D]>2, x = ConstantArray[0,{k,2}], x = ConstantArray[0,{k}]];
	Do[
		time = 0;
	    (* play state transitions *)
        While[state<=NN,
            time -= Log[RandomReal[]] sojourn[[state]];
            r = RandomReal[];
            nstate = 1;
            While[nextpr[[state,nstate]]<=r, nstate++];
            state = nstate;
        ];
		If[Length[D]>2, 
			x[[n,1]]=time; 
			x[[n,2]]=Ceiling[state/NN]-1;
		,
			x[[n]] = time;
		];
	,{n,k}];
	Return[x];
];


End[(* Private *)];
EndPackage[];
