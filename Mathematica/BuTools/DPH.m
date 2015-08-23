(* ::Package:: *)

(*
   BuTools DPH Package
*)

BeginPackage["BuTools`DPH`"];
AcyclicDPHFromMG::usage = "{beta, B} = AcyclicDPHFromMG[alpha, A, precision]: Transforms a matrix-geometric representation to an acyclic DPH representation of the same size, if possible.";
CanonicalFromDPH2::usage = "{beta, B} = CanonicalFromDPH2[alpha, A, prec]: Returns the canonical form of an order-2 discrete phase-type distribution.";
CanonicalFromDPH3::usage = "{beta, B} = CanonicalFromDPH3[alpha, A, prec]: Returns the canonical form of an order-3 discrete phase-type distribution.";
CdfFromMG::usage = "cdf = CdfFromMG[alpha, A, x, prec]: Returns the cummulative distribution function of a matrix-geometric distribution.";
CdfFromDPH::usage = "cdf = CdfFromDPH[alpha, A, x, prec]: Returns the cummulative distribution function of a discrete phase-type distribution.";
CheckMGRepresentation::usage = "r = CheckMGRepresentation[alpha, A, prec]: Checks if the given vector and matrix define a valid matrix-geometric representation.";
CheckDPHRepresentation::usage = "r = CheckDPHRepresentation[alpha, A, prec]: Checks if the given vector and matrix define a valid discrete phase-type representation.";
DPH2From3Moments::usage = "{alpha, A} = DPH2From3Moments[moms, prec]: Returns an order-2 discrete phase-type distribution which has the same 3 moments as given.";
DPH3From5Moments::usage = "{alpha, A} = DPH3From5Moments[moms, prec]: Returns an order-3 discrete phase-type distribution which has the same 5 moments as given.";
DPHFromMG::usage = "{beta, B} = DPHFromMG[alpha, A, precision]: Obtains a Markovian representation of a matrix geometric distribution of the same size, if possible.";
MGFromMoments::usage = "{alpha, A} = MGFromMoments[moms]: Creates a matrix-geometric distribution that has the same moments as given.";
MGOrderFromMoments::usage = "alma";
MomentsFromMG::usage = "moms = MomentsFromMG[alpha, A, K, prec]: Returns the moments of a matrix geometric distribution.";
MomentsFromDPH::usage = "moms = MomentsFromDPH[alpha, A, K, prec]: Returns the moments of a discrete phase-type distribution.";
PmfFromMG::usage = "pmf = PmfFromMG[alpha, A, x, prec]: Returns the probability mass function of a matrix-geometric distribution.";
PmfFromDPH::usage = "pmf = PmfFromDPH[alpha, A, x, prec]: Returns the probability mass function of a discrete phase-type distribution.";
RandomDPH::usage = "{alpha, A} = RandomDPH[order, mean, zeroEntries, maxTrials, prec]: Returns a random discrete phase-type distribution with a given mean value.";
SamplesFromDPH::usage = "x = SamplesFromDPH[alpha, A, K, prec]: Generates random samples from a discrete phase-type distribution.";


Begin["`Private`"];


Needs["BuTools`RepTrans`"];
Needs["BuTools`Moments`"];
Needs["BuTools`MC`"];
Needs["BuTools`PH`"];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckInput"]],BuTools`CheckInput=True];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckPrecision"]],BuTools`CheckPrecision=N[10^-12]];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`Verbose"]],BuTools`Verbose=False];


AcyclicDPHFromMG[alpha_, A_, precision_:N[10^-14]]:=
Module[{G,T,gamma,beta,B,lambda,lambda2,mx},
	If[BuTools`CheckInput && Not[CheckMGRepresentation[alpha, A]],Throw["AcyclicDPHFromMG: Input is not a valid MG distribution!"]];
	lambda = Eigenvalues[A];
    If[Max[Abs[Im[lambda]]]>precision,Throw["AcyclicDPHFromMG: The input matrix has complex eigenvalue!"]];
	lambda2 = Sort[Eigenvalues[A],Abs[#1]>Abs[#2]||(Abs[#1]==Abs[#2] && Re[#1]>Re[#2])||(Abs[#1]==Abs[#2]&&Re[#1]==Re[#2]&&Im[#1]>=Im[#2])&];    
    mx = DiagonalMatrix[lambda2]+DiagonalMatrix[1-lambda2[[;;-2]], 1];
    T = SimilarityMatrix[A, mx];
    beta = alpha.T;
    B = mx;
	If[Not[CheckDPHRepresentation[beta, B, precision]], Throw["AcyclicDPHFromMG: No acyclic representation found!"]];
	Return[{beta,B}];
];


CanonicalFromDPH2[alpha_, A_, prec_:N[10^-14]]:=
Module[{\[Lambda]1,\[Lambda]2,e,p1,\[Alpha]v,Av,\[Delta]1,\[Delta]2},
    If[BuTools`CheckInput && Not[CheckMGRepresentation[alpha,A]], Throw["CanonicalFromDPH2: Input is not a valid MG distribution!"]];
    If[Dimensions[A][[1]]!=2, Throw["CanonicalFromDPH2: Dimension is not 2!"]];
    
	{\[Lambda]1,\[Lambda]2}=Sort[Eigenvalues[A],Abs[#1]>Abs[#2]||(Abs[#1]==Abs[#2] && Re[#1]>Re[#2])||(Abs[#1]==Abs[#2]&&Re[#1]==Re[#2]&&Im[#1]>=Im[#2])&];
	e={{1},{1}};
	p1=(alpha.(e-A.e))[[1]];

	Which[\[Lambda]1>0 && \[Lambda]2>=0 && \[Lambda]1!=\[Lambda]2,
		\[Delta]1=(1-\[Lambda]1)(1-p1-\[Lambda]2)/(\[Lambda]1-\[Lambda]2);
		\[Delta]2=p1-\[Delta]1;
		\[Alpha]v={\[Delta]1 (\[Lambda]1-\[Lambda]2)/((1-\[Lambda]1)(1-\[Lambda]2)),(\[Delta]1+\[Delta]2)/(1-\[Lambda]2)};
		Av={{\[Lambda]1,1-\[Lambda]1},{0,\[Lambda]2}};
	,
		\[Lambda]1>0 && \[Lambda]1==\[Lambda]2,(*Same eigenvalues*)
		\[Delta]2=p1;
		\[Delta]1=(1-\[Lambda]1)(1-\[Delta]2-\[Lambda]1)/\[Lambda]1;
		\[Alpha]v={\[Delta]1 \[Lambda]1/(1-\[Lambda]1)^2,\[Delta]2/(1-\[Lambda]1)};
		Av={{\[Lambda]1,1-\[Lambda]1},{0,\[Lambda]1}};
	,
		\[Lambda]1>0,(*Negative eigenvalues*)
		\[Delta]1=(1-\[Lambda]1)(1-p1-\[Lambda]2)/(\[Lambda]1-\[Lambda]2);
		\[Delta]2=p1-\[Delta]1;
		\[Alpha]v={(\[Delta]1 \[Lambda]1+\[Delta]2 \[Lambda]2)/((1-\[Lambda]1)(1-\[Lambda]2)),(\[Delta]1+\[Delta]2)(1-\[Lambda]1-\[Lambda]2)/((1-\[Lambda]1)(1-\[Lambda]2))};
		Av={{\[Lambda]1+\[Lambda]2,1-\[Lambda]1-\[Lambda]2},{\[Lambda]1 \[Lambda]2/(\[Lambda]1+\[Lambda]2-1),0}};
	];
	Return[{\[Alpha]v,Av}];
];


CanonicalFromDPH3[alpha_, A_, prec_:N[10^-14]]:=
Module[{\[Lambda]1,\[Lambda]2,\[Lambda]3,\[Alpha]v,Av,A2,\[Delta]1,\[Delta]2,\[Delta]3,b1,b2,b3,B,e,a0,a1,a2,x,x1,x2,x3,p1,p2,
M,m1,m2,m3,m33,res,i},
    If[BuTools`CheckInput && Not[CheckMGRepresentation[alpha,A]], Throw["CanonicalFromDPH3: Input is not a valid MG distribution!"]];
    If[Dimensions[A][[1]]!=3, Throw["CanonicalFromDPH3: Dimension is not 3!"]];
    
	{\[Lambda]1,\[Lambda]2,\[Lambda]3}=Sort[Eigenvalues[A],Abs[#1]>Abs[#2]||(Abs[#1]==Abs[#2] && Re[#1]>Re[#2])||(Abs[#1]==Abs[#2]&&Re[#1]==Re[#2]&&Im[#1]>=Im[#2])&];
	If[BuTools`Verbose,Print["Ordered eigenvalues:",{\[Lambda]1,\[Lambda]2,\[Lambda]3}]];
	a0=-\[Lambda]1 \[Lambda]2 \[Lambda]3;
	a1=\[Lambda]1 \[Lambda]2+\[Lambda]1 \[Lambda]3+\[Lambda]2 \[Lambda]3;
	a2=-\[Lambda]1-\[Lambda]2-\[Lambda]3;
	e=Table[1,{3}];

	Which[Re[\[Lambda]1]<0,
		Throw["DPH3Canonical: The largest eigenvalue has negative real part!"],
		Not[Or @@ Negative[Re[{\[Lambda]1,\[Lambda]2,\[Lambda]3}]]],
		(* PPP *) 
		{\[Alpha]v,A2}=CanonicalFromPH3[alpha,A-IdentityMatrix[3],prec];
		Av=A2+IdentityMatrix[3];,
		And[Re[\[Lambda]2]>=0,Re[\[Lambda]3]<0],
		(* PPN *)
		x1=\[Lambda]1;
		x2=\[Lambda]2+\[Lambda]3;
		x3=\[Lambda]2 \[Lambda]3/(\[Lambda]2+\[Lambda]3-1);
		Av={{x1,1-x1,0},{0,x2,1-x2},{0,x3,0}};
		b3=1/(1-x3)(e-A.e);
		b2=1/(1-x2)A.b3;
		b1=e-b2-b3;
		B=Transpose[{b1,b2,b3}];
		\[Alpha]v=alpha.B;,
		And[Re[\[Lambda]3]>=prec,Re[\[Lambda]2]<0],
		(* PNP+ *)
		x1=\[Lambda]3;
		x2=\[Lambda]1+\[Lambda]2;
		x3=\[Lambda]1 \[Lambda]2/(\[Lambda]1+\[Lambda]2-1);
		p1=alpha.(e-A.e);
		p2=alpha.A.(e-A.e);
		\[Delta]1=(1-\[Lambda]1)((1-\[Lambda]2)(1-\[Lambda]3)+(-1+\[Lambda]2+\[Lambda]3)p1-p2)/((\[Lambda]1-\[Lambda]2)(\[Lambda]1-\[Lambda]3));
		\[Delta]2=(\[Lambda]2-1)((1-\[Lambda]1)(1-\[Lambda]3)+(-1+\[Lambda]1+\[Lambda]3)p1-p2)/((\[Lambda]1-\[Lambda]2)(\[Lambda]2-\[Lambda]3));
		\[Delta]3=(\[Lambda]3-1)((1-\[Lambda]1)(1-\[Lambda]2)+(-1+\[Lambda]1+\[Lambda]2)p1-p2)/((\[Lambda]2-\[Lambda]3)(\[Lambda]3-\[Lambda]1));
		Av={{x1,0,0},{0,x2,1-x2},{0,x3,0}};
		\[Alpha]v={\[Delta]3/(1-\[Lambda]3),(\[Delta]1 \[Lambda]1+\[Delta]2 \[Lambda]2)/((1-\[Lambda]1)(1-\[Lambda]2)),(\[Delta]1+\[Delta]2)(1-\[Lambda]1-\[Lambda]2)/((1-\[Lambda]1)(1-\[Lambda]2))};
		If[Min[alpha,A]<-prec,
			(* PNP *)
			x1=-a2;
			x2=(a0-a1 a2)/(a2 (1+a2));
			x3=a0 (1+a2)/(a0-a2-a1 a2-a2^2);
			Av={{x1,1-x1,0},{x2,0,1-x2},{0,x3,0}};
			b3=1/(1-x3)(e-A.e);
			b2=1/(1-x2)A.b3;
			b1=e-b2-b3;
			B=Transpose[{b1,b2,b3}];
			\[Alpha]v=alpha.B;
			If[Min[alpha.B]<-prec,
				(* Set the initial vector first element > 0 *)
				M=({{m1, 1-m1, 0},{m2, 0, 1-m2},{0, m3, m33}});
			    res=Solve[Eigenvalues[M]=={\[Lambda]2,\[Lambda]3,\[Lambda]1},{m1,m2,m3}];
				{m1,m2,m3}={m1,m2,m3}/.res[[1]];
				b3=(IdentityMatrix[3]-A).e/(1-m3-m33);
				b2=(- m33 IdentityMatrix[3]+A).b3/(1-m2);
				b1=(- m3 b3+A.b2)/(1-m1);
				B=Transpose[{b1,b2,b3}];
				\[Alpha]v=alpha.B;
				x=0;
				While[x<=1,
					If[\[Alpha]v[[1]]/.m33->x >= 0 && m1/.m33->x >= 0 && m2/.m33->x >= 0 && m3/.m33->x >= 0,Break[]];
					x=x+1/100;
				];
				m33=x;
				Av=M;
				If[Min[\[Alpha]v,Av]<0,Print["CanonicalFromDPH3: Unhandled PNP case! Input:",alpha,A];Throw["CanonicalFromDPH3: PNP ERROR!"]];
			]
		],
		And[Re[\[Lambda]2]<0,Re[\[Lambda]3]<0],
		(* PNN *)
		If[!Element[{\[Lambda]1,\[Lambda]2,\[Lambda]3},Reals]&&Abs[\[Lambda]2]^2>2 \[Lambda]1 (-Re[\[Lambda]2]),
			{\[Alpha]v,A2}=PH3Canonical[alpha,A-IdentityMatrix[3],prec];
			Av=A2+IdentityMatrix[3];
		,
			x1=-a2;
			x2=-a1/(1+a2);
			x3=-a0/(1+a1+a2);
			Av={{x1,1-x1,0},{x2,0,1-x2},{x3,0,0}};
			b3=1/(1-x3)(e-A.e);
			b2=1/(1-x2)A.b3;
			b1=e-b2-b3;
			B=Transpose[{b1,b2,b3}];
			\[Alpha]v=alpha.B;
			If[Max[Abs[Im[\[Alpha]v]]]>prec,Print["Input:",alpha,","A,"Out vector:",\[Alpha]v];Throw["CanonicalFromDPH3: canonical vector is complex! Input:"],\[Alpha]v=Re[\[Alpha]v]];
			If[Max[Abs[Im[Av]]]>prec,Print["Input:",alpha,","A,"Out matrix:",Av];Throw["CanonicalFromDPH3: canonical matrix is complex!"],Av=Re[Av]];
		];
		,
		True,
		Print["Unhandled"];
	];
	Return[{\[Alpha]v,Av}];
];


CdfFromMG[alpha_, A_, x_]:=(
    If[BuTools`CheckInput && Not[CheckMGRepresentation[alpha, A]],Throw["CdfFromMG: Input is not a valid MG distribution!"]];
	If[ListQ[x],
		Return[Table[1 - Total[If[xv==0,alpha,alpha.MatrixPower[A,xv]]],{xv,x}]]
	, Return[1 - Total[If[x==0,alpha,alpha.MatrixPower[A,x]]]]]);


CdfFromDPH[alpha_, A_, x_]:=(
    If[BuTools`CheckInput && Not[CheckDPHRepresentation[alpha, A]],Throw["CdfFromDPH: Input is not a valid DPH distribution!"]];
	Return[CdfFromMG[alpha,A,x]]);


CheckDPHRepresentation[alpha_, A_, precision_:Null]:=
Module[{prec},
	If[Not[NumericQ[precision]], prec=BuTools`CheckPrecision, prec=precision];
	If[Dimensions[alpha][[1]]!=Dimensions[A][[1]],
		If[BuTools`Verbose,Print["CheckDPHRepresentation: alpha and A have different sizes!"]];
		Return[False]
	];
	Return[CheckProbVector[alpha,True,prec] && CheckProbMatrix[A,True,prec]];
];


CheckMGRepresentation[alpha_, A_,precision_ :Null]:=
Module[{eig,prec},

	If[Not[NumericQ[precision]], prec=BuTools`CheckPrecision, prec=precision];
	If[ Dimensions[A][[1]]!=Dimensions[A][[2]],
		If[BuTools`Verbose==True,Print["CheckMGRepresentation: The matrix is not a quadratic matrix!"]];
		Return[False]
	];

	If[ Dimensions[alpha][[1]]!=Dimensions[A][[1]],
		If[BuTools`Verbose==True,Print["CheckMGRepresentation: The vector and the matrix have different sizes!"]];
		Return[False]
	];

	If[Total[alpha]<-prec || Total[alpha]>1+prec,
		If[BuTools`Verbose==True,Print["CheckMGRepresentation: The sum of the vector is not 1 (precision:",prec,")!"]];
		Return[False]
	];

	eig=Sort[Eigenvalues[A],Abs[#1]>Abs[#2]||(Abs[#1]==Abs[#2] && Re[#1]>Re[#2])||(Abs[#1]==Abs[#2]&&Re[#1]==Re[#2]&&Im[#1]>=Im[#2])&];

	If[Abs[Im[eig[[1]]]] > prec,
		If[BuTools`Verbose==True,Print["CheckMGRepresentation: The largest eigenvalue of the matrix is complex!"]];
        Return[False]
	];

	If[eig[[1]]> 1+prec,
		If[BuTools`Verbose==True,Print["CheckMGRepresentation: The largest eigenvalue of the matrix is greater than 1 (precision:",prec,")!"]];
		Return[False]
	];

	If[Abs[eig[[1]]]==Abs[eig[[2]]],
		If[BuTools`Verbose==True,Print["CheckMGRepresentation Warning: There are more than one eigenvalue with the same absolute value as the largest eigenvalue!"]]
	];
	Return[True];
];


DPH2From3Moments[moms_]:=
Module[{beta,B},
    {beta, B} = MGFromMoments[moms[[1;;3]]];
    Return[CanonicalFromDPH2[beta,B]];
];


DPH3From5Moments[moms_, prec_:N[10^-14]]:=
Module[{beta,B},
    {beta, B} = MGFromMoments[moms[[1;;5]]];
    Return[CanonicalFromDPH3[beta,B,prec]];
];


DPHFromMG[alpha_, A_, precision_:N[10^-14]]:=
Module[{Transfun,Evalfun,nrep},
    Transfun[orep_, B_]:=
        Return[{orep[[1]].B, Inverse[B].orep[[2]].B}];
  
    Evalfun[orep_, k_:0]:=Module[{ao,Ao,av,Ad},
        ao = orep[[1]];
        Ao = orep[[2]];
        av = 1-Total[Ao,{2}];
        Ad = Ao - DiagonalMatrix[Diagonal[Ao]];
        If[Mod[k,2] == 0,
            Return[-Min[ao, av, Ad]];
        ,
            Return[-Total[Select[ao,Negative]] - Total[Select[av,Negative]] - Total[Select[Flatten[Ad],Negative]]];
        ]];
   
    If[BuTools`CheckInput && Not[CheckMGRepresentation[alpha, A]],Throw["DPHFromMG: Input is not a valid MG distribution!"]];
    nrep = FindMarkovianRepresentation[{alpha, A}, Transfun, Evalfun, precision];
    Return[nrep];
];


MGFromMoments[moms_]:=
Module[{rfmom,\[Alpha],A,C,vlist},
	rfmom=ReducedMomsFromMoms[FactorialMomsFromMoms[moms]];
	rfmom=Join[{1},rfmom];
	vlist=Table[k!Sum[(-1)^i rfmom[[i+1]],{i,0,k}],{k,1,Length[moms]}];
	{\[Alpha],C}=MEFromMoments[vlist];
	A=Inverse[C].Inverse[Inverse[C]+IdentityMatrix[Dimensions[C][[1]]]];
	Return[{\[Alpha],A}];
];


MGOrderFromMoments[moms_,prec_ : N[10^-10]]:=
Module[{size,rmoms,n,hankel,order},
	
	size = Floor[(Length[moms]+1)/2];
	rmoms = Prepend[ReducedMomsFromMoms[FactorialMomsFromMoms[moms]],1];
	order = size;
	For[n=1,n<=size,n++,
		hankel=Table[rmoms[[i+j-1]],{i,1,n},{j,1,n}];
		If[Abs[Det[hankel]] < prec, order=n-1; Break[]];
	];
	Return[order];
];


MomentsFromDPH[alpha_,A_,K_:0]:=(
    If[BuTools`CheckInput && Not[CheckDPHRepresentation[alpha, A]],Throw["MomentsFromDPH: Input is not a valid DPH representation!"]];
    MomentsFromMG[alpha,A,K]);


MomentsFromMG[alpha_,A_,K_:0]:=
Module[{KK,fmoms,iA},
    If[BuTools`CheckInput && Not[CheckMGRepresentation[alpha, A]],Throw["MomentsFromMG: Input is not a valid MG representation!"]];
    If[K==0,KK=2 Length[alpha]-1,KK=K];
	iA=Inverse[IdentityMatrix[Dimensions[A][[1]]]-A];
	fmoms=Table[i! Total[alpha.MatrixPower[iA,i].If[i==1,IdentityMatrix[Length[A]],MatrixPower[A,i-1]]],{i,1,KK}];
	Return[MomsFromFactorialMoms[fmoms]];
]


PmfFromMG[alpha_, A_, x_]:=
Module[{a},
    If[BuTools`CheckInput && Not[CheckMGRepresentation[alpha, A]],Throw["PmfFromMG: Input is not a valid MG distribution!"]];
	a=1-Total[A,{2}];
	If[ListQ[x],
		Return[Table[If[xv==0,1-Total[alpha],If[xv==1,alpha,alpha.MatrixPower[A, xv-1]].a],{xv,x}]]
	, Return[If[x==0,1-Total[alpha],If[x==1,alpha,alpha.MatrixPower[A, x-1]].a]];
	];
];


PmfFromDPH[alpha_, A_, x_]:=(
    If[BuTools`CheckInput && Not[CheckDPHRepresentation[alpha, A]],Throw["PmfFromDPH: Input is not a valid DPH distribution!"]];
	Return[PmfFromMG[alpha,A,x]]);


RandomDPH[order_,mean_:10,zeroEntries_:0,maxTrials_:1000,prec_: N[10^-7]]:=
Module[{trials,zeros,numZeros,idx,zeroInRow,Row,B,sumRow,elem,vector,actualZeros,d,m},
If[zeroEntries>(order+1)(order-1),Throw["RandomDPH: You have given too many zeros! Try to decrease the zero entries number!"]];
actualZeros=zeroEntries;
While[actualZeros>=0,
trials=1;
While[trials < maxTrials,
zeros=Table[0,{i,order+1}];
numZeros=0;
idx=1;
While[numZeros<actualZeros && idx<order+1,
zeroInRow=RandomInteger[Min[actualZeros-numZeros,order-1]];
zeros[[idx]]=zeroInRow;
numZeros=numZeros+zeroInRow;
idx++;
];
zeros[[order+1]]=actualZeros-numZeros;
If[zeros[[order]]>=order,trials++;Continue[]];(*Too many zeros in the last row*)
B=Table[0,{i,order},{j,order}];
For[idx=1,idx<=order,idx++,
sumRow=0;
Row=Table[elem=Random[Real,1-sumRow];sumRow=sumRow+elem;elem,{i,order-zeros[[idx]]}];
Row=Join[Row,{1-sumRow}];
Row=Join[Row,Table[0,{i,zeros[[idx]]}]];
Row=Row[[RandomSample[Range[order+1],order+1]]];
B[[idx]]=Row[[1;;order]];
];
If[MatrixRank[IdentityMatrix[order]-B]!=order,trials++;Continue[]];
sumRow=0;
vector=Table[elem=Random[Real,1-sumRow];sumRow=sumRow+elem;elem,{i,order-zeros[[order+1]]-1}];
vector=Join[vector,{1-sumRow}];
vector=Join[vector,Table[0,{i,zeros[[order+1]]}]];
vector=vector[[RandomSample[Range[order],order]]];
If[Min[Abs[vector.Inverse[IdentityMatrix[order]-B]]] < prec,trials++;Continue[]];
d = Table[Random[Real],{order}];
(* scale to the mean value *)
m = MomentsFromDPH[vector, DiagonalMatrix[1-d].B+DiagonalMatrix[d], 1];
d = 1 - (1-d) m[[1]]/mean;
B = DiagonalMatrix[1-d].B+DiagonalMatrix[d];
If[Not[CheckDPHRepresentation[vector,B,prec]],Continue[]];
If[zeroEntries>actualZeros,
Print["RandomDPH: Number of zero entries are different! Given:",zeroEntries,", in the returned representation:",actualZeros]];
Return[{vector,B}];
];
--actualZeros;
];
];


SamplesFromDPH[a_,A_,k_]:=
Module[{NN,cummInitial,sojourn, nextpr,logp,genSamples},
	If[BuTools`CheckInput && Not[CheckDPHRepresentation[a,A]],Throw["SamplesFromDPH: input is not a valid DPH representation!"]];
    (* auxilary variables*)
    NN = Length[a];
    cummInitial = Accumulate[a];
	logp = Log[Diagonal[A]];
    sojourn = 1/(1-Diagonal[A]);
    nextpr = DiagonalMatrix[sojourn].A;
    nextpr = nextpr - DiagonalMatrix[Diagonal[nextpr]];
    nextpr = Join[nextpr, Transpose[{1-Total[nextpr,{2}]}],2];
    nextpr = Transpose[Accumulate[Transpose[nextpr]]];
    
    genSamples=Compile[{{NN,_Integer},{logp,_Real,1},{cummInitial,_Real,1},{sojourn,_Real,1},{nextpr,_Real,2}},
	Module[{time,rr,state,nstate,x},
		x = Table[0,{k}];
		Do[
			time = 0;
			(* draw initial distribution *)
			rr = RandomReal[];
			state = 1;
			While[cummInitial[[state]]<=rr, state++];
			(* play state transitions *)
			While[state<=NN,
				time = time + 1 + Floor[Log[RandomReal[]] / logp[[state]]];
				rr = RandomReal[];
				nstate = 1;
				While[nextpr[[state,nstate]]<=rr, nstate++];
				state = nstate;
			];
			x[[n]] = time;
		,{n,k}];
		Return[x]
	]];
	Return[genSamples[NN,logp,cummInitial,sojourn,nextpr]]
];


End[(*Private*)];
EndPackage[];
