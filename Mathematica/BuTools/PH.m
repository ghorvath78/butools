(* ::Package:: *)

(*
   BuTools_PH Package
*)

BeginPackage["BuTools`PH`"];
AcyclicPHFromME::usage = "{beta, B} = AcyclicPHFromME[alpha, A, maxSize, precision]: Transforms an arbitrary matrix-exponential representation to an acyclic phase-type representation. ";
APH2ndMomentLowerBound::usage = "m2 = APH2ndMomentLowerBound[m1, n]: Returns the lower bound of the second moment of acyclic phase-type (APH) distributions of order n.";
APH3rdMomentLowerBound::usage = "m3 = APH3rdMomentLowerBound[m1, m2, n]: Returns the lower bound of the third moment of acyclic phase-type (APH) distributions of order n.";
APH3rdMomentUpperBound::usage = "m3 = APH3rdMomentUpperBound[m1, m2, n]: Returns the upper bound of the third moment of acyclic phase-type (APH) distributions of order n.";
APHFrom2Moments::usage = ": Returns an acyclic PH which has the same 2 moments as given.";
APHFrom3Moments::usage = "{alpha, A} = APHFrom3Moments[moms, maxSize]: Returns an acyclic PH which has the same 3 moments as given.";
CanonicalFromPH2::usage = "{beta, B} = CanonicalFromPH2[alpha, A, prec]: Returns the canonical form of an order-2 phase-type distribution.";
CanonicalFromPH3::usage = "{beta, B} = CanonicalFromPH3[alpha, A, prec]: Returns the canonical form of an order-3 phase-type distribution.";
CdfFromME::usage = "cdf = CdfFromME[alpha, A, x]: Returns the cummulative distribution function of a matrix-exponential distribution.";
CdfFromPH::usage = "cdf = CdfFromPH[alpha, A, x]: Returns the cummulative distribution function of a continuous phase-type distribution.";
CheckMEPositiveDensity::usage = "r = CheckMEPositiveDensity[alpha, A, maxSize, prec]: Checks if the given matrix-exponential distribution has positive density.";
CheckMERepresentation::usage = "r = CheckMERepresentation[alpha, A, prec]: Checks if the given vector and matrix define a valid matrix-exponential representation.";
CheckPHRepresentation::usage = "r = CheckPHRepresentation[alpha, A, prec]: Checks if the given vector and matrix define a valid phase-type representation.";
IntervalPdfFromPH::usage = "{x,y} = IntervalPdfFromPH[alpha, A, intBounds, prec]: Returns the approximate probability density function of a continuous phase-type distribution, based on the probability of falling into intervals.";
MEFromMoments::usage = "{alpha, A} = MEFromMoments[moms]: Creates a matrix-exponential distribution that has the same moments as given.     ";
MEOrder::usage = "order = MEOrder[alpha, A, kind, prec]: Returns the order of the ME distribution (which is not necessarily equal to the size of the representation).";
MEOrderFromMoments::usage = "order = MEOrderFromMoments[moms, prec]: Returns the order of ME distribution that can realize the given moments.";
MinimalRepFromME::usage = "{beta, B} = MinimalRepFromME[alpha, A, how, precision]: Returns the minimal representation of the given ME distribution.";
MomentsFromME::usage = "moms = MomentsFromME[alpha, A, K, prec]: Returns the moments of a matrix-exponential distribution.";
MomentsFromPH::usage = "moms = MomentsFromPH[alpha, A, K, prec]: Returns the moments of a continuous phase-type distribution.";
MonocyclicPHFromME::usage = "{beta, B} = MonocyclicPHFromME[alpha, A, maxSize, precision]: Transforms an arbitrary matrix-exponential representation to a Markovian monocyclic representation.";
PdfFromME::usage = "pdf = PdfFromME[alpha, A, x, prec]: Returns the probability density function of a matrix-exponential distribution.";
PdfFromPH::usage = "pdf = PdfFromPH[alpha, A, x, prec]: Returns the probability density function of a continuous phase-type distribution.";
PH2From3Moments::usage = "{alpha, A} = PH2From3Moments[moms, prec]: Returns a PH(2) which has the same 3 moments as given.";
PH3From5Moments::usage = "{alpha, A} = PH3From5Moments[moms]: Returns a PH(3) which has the same 5 moments as given.";
PHFromME::usage = "{beta, B} = PHFromME[alpha, A, precision]: Obtains a Markovian representation of a matrix exponential distribution of the same size, if possible.";
RandomPH::usage = "{alpha, A} = RandomPH[order, mean, zeroEntries, maxTrials, prec]: Returns a random phase-type distribution with a given order.";
SamplesFromPH::usage = "x = SamplesFromPH[alpha, A, K, prec]: Generates random samples from a phase-type distribution.";
TestPHPackage::usage = "TestPHPackage[] : Executes various tests to check the functions of the PH package.";


Begin["`Private`"];


Needs["BuTools`RepTrans`"];
Needs["BuTools`Moments`"];
Needs["BuTools`MC`"];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckInput"]],BuTools`CheckInput=True];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckPrecision"]],BuTools`CheckPrecision=N[10^-12]];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`Verbose"]],BuTools`Verbose=False];


AcyclicPHFromME[alpha_, A_, maxSize_:100, precision_:N[10^-14]]:=
Module[{G,T,gamma,beta,B},
	If[BuTools`CheckInput && Not[CheckMERepresentation[alpha, A]],Throw["AcyclicPHFromME: Input is not a valid ME distribution!"]];
	G = TransformToAcyclic[A, maxSize, precision];
    T = SimilarityMatrix[A, G];
    gamma = Re[alpha.T];
    If[Min[gamma]>-precision; 
		beta=gamma;B = G,
        {beta, B} = ExtendToMarkovian[gamma, G, maxSize, precision];
	];
    If[Not[CheckPHRepresentation[beta, B, precision]],Throw["AcyclicPHFromME: No acyclic representation found up to the given size and precision!"]];
	Return[{beta,B}];
];


APH2ndMomentLowerBound[m1_,n_]:= m1^2 (n+1)/n;


APH3rdMomentLowerBound[m1_,m2_,n_]:=
Module[{n2,p,l,a},
    n2 = m2 / m1 / m1;
    If[n2<(n+1)/n,
        Return[Infinity]
    , If[n2<(n+4)/(n+1),
        p = ((n+1)*(n2-2)) / (3*n2*(n-1)) * ((-2*Sqrt[n+1]) / Sqrt[-3*n*n2+4*n+4] -1);
        a = (n2-2) / (p*(1-n2) + Sqrt[p*p+p*n*(n2-2)/(n-1)]);
        l = ((3+a)*(n-1)+2*a) / ((n-1)*(1+a*p)) - (2*a*(n+1)) / (2*(n-1)+a*p*(n*a+2*n-2));
        Return[Re[l] * m1 * m2]
	  , Return[(n+1.0)/n * n2 * m1 * m2]
	  ]
	];
];


APH3rdMomentUpperBound[m1_,m2_,n_]:=
Module[{n2},
    n2 = m2 / m1 / m1;
    If[n2<(n+1)/n,
        Return[-Infinity];
    , If[n2<=n/(n-1),
        Return[m1 * m2 * (2*(n-2)*(n*n2-n-1)*Sqrt[1+(n*(n2-2))/(n-1)] + (n+2)*(3*n*n2-2*n-2)) / (n*n*n2)];
	,   Return[Infinity];
	  ];
	];
];


APHFrom2Moments[moms_]:=
Module[{cv2,lambda,p,n,A,alpha},
    cv2 = moms[[2]]/moms[[1]]^2 - 1;
    lambda = 1 / moms[[1]];
    n = Max[Ceiling[1/cv2], 2];
    p = 1 / (cv2 + 1 + (cv2-1)/(n-1));
	A=ConstantArray[0,{n,n}];
    A = -lambda*p*n * IdentityMatrix[n];
    Do[A[[i,i+1]] = -A[[i,i]],{i,n-1}];
    A[[n,n]] = -lambda*n;
    alpha = ConstantArray[0,{n}];
    alpha[[1]] = p;
    alpha[[n]] = 1 - p;
	Return[{alpha,A}];
];


APHFrom3Moments[moms_, maxSize_:100, prec_:N[10^-14]]:=
Module[{m1,m2,m3,n,n1,n2,n3,b,a,p,lambda,mu,c0,c1,c2,c3,c4,fs,alpha,A},
	{m1,m2,m3} = moms;
    
    (* detect number of phases needed *)
    n = 2;
    While[n<maxSize && (APH2ndMomentLowerBound[m1, n] > m2 || APH3rdMomentLowerBound[m1, m2, n] >= m3 || APH3rdMomentUpperBound[m1, m2, n] <= m3), n++];
        
    (* if PH is too large, adjust moment to bounds *)
    If[APH2ndMomentLowerBound[m1, n] > m2, 
        m2 = APH2ndMomentLowerBound[m1, n];
		Print["APHFrom3Moments: moments are not feasible with maxSize phases. Adjusting second moment."];
    ];
    If[APH3rdMomentLowerBound[m1, m2, n] > m3,
        m3 = APH3rdMomentLowerBound[m1, m2, n];
		Print["APHFrom3Moments: moments are not feasible with maxSize phases. Adjusting third moment."];
	];
    If[APH3rdMomentUpperBound[m1, m2, n] < m3,
        m3 = APH3rdMomentUpperBound[m1, m2, n];
		Print["APHFrom3Moments: moments are not feasible with maxSize phases. Adjusting third moment."];
    ];
	(* compute normalized moments *)
    {n1,n2,n3} = NormMomsFromMoms[{m1,m2,m3}];
    If[n2>2 || n3 < 2n2 - 1,
        b = Re[2*(4-n*(3*n2-4)) / (n2*(4+n-n*n3) + Sqrt[n*n2]*Sqrt[12*n2^2*(n+1)+16*n3*(n+1)+n2*(n*(n3-15)*(n3+1)-8*(n3+3))])];
        a = (b*n2-2)*(n-1)*b / (b-1) / n;
        p = (b-1) / a;
        lambda = (p*a+1) / n1;
        mu = (n-1)*lambda / a;
        (* construct representation *)
        alpha = ConstantArray[0,{n}];
        alpha[[1]] = p;
        alpha[[n]] = 1.0-p;
        A = ConstantArray[0,{n,n}];
        A[[n,n]] = -lambda;
		Do[A[[i,i]]=-mu;A[[i,i+1]]=mu,{i,1,n-1}];
        Return[{alpha,A}];
    ,
        c4 = n2*(3*n2-2*n3)*(n-1)^2;
        c3 = 2*n2*(n3-3)*(n-1)^2;
        c2 = 6*(n-1)*(n-n2);
        c1 = 4*n*(2-n);
        c0 = n*(n-2);
		fs=Map[#1[[1,2]]&,Solve[c4 f^4+c3 f^3+c2 f^2+c1 f+c0==0,f]];
		Do[
            If[Abs[(n-1)*(n2*f^2-2*f+2)-n]<prec, Continue[]];
            a = 2*(f-1)*(n-1) / ((n-1)*(n2*f^2-2*f+2)-n);
            p = (f-1)*a;
            lambda = (a+p) / n1;
            mu = (n-1) / (n1 - p/lambda);
            If[ Element[p,Reals] && Element[lambda,Reals] && Element[mu,Reals] && p>=0 && p<=1 && lambda>0 && mu>0,	
				alpha = ConstantArray[0,{n}];
                alpha[[1]] = p;
                alpha[[2]] = 1-p;
                A = ConstantArray[0,{n,n}];
                A[[1,1]] = -lambda;
                A[[1,2]] = lambda;
				Do[A[[j,j]]=-mu; If[j<n,A[[j,j+1]]=mu],{j,2,n}];
                Return[{alpha,A}];
            ];
		,{f,fs}];
		If[Not[ListQ[alpha]],Throw["No APH found for the given 3 moments!"], Return[{alpha,A}]];
   ];
];


CanonicalFromPH2[alpha_, A_, prec_:N[10^-14]]:=(
    If[BuTools`CheckInput && Not[CheckMERepresentation[alpha,A]], Throw["CanonicalFromPH2: Input is not a valid ME distribution!"]];
    If[Dimensions[A][[1]]!=2, Throw["CanonicalFromPH2: Dimension is not 2!"]];
    PH2From3Moments[MomentsFromME[alpha, A, 3], prec]);


CanonicalFromPH3[alpha_, A_, prec_:N[10^-14]]:=(
    If[BuTools`CheckInput && Not[CheckMERepresentation[alpha,A]], Throw["CanonicalFromPH3: Input is not a valid ME distribution!"]];
    If[Dimensions[A][[1]]!=3, Throw["CanonicalFromPH3: Dimension is not 3!"]];
    PH3From5Moments[MomentsFromME[alpha, A, 5], prec]);


CdfFromME[alpha_, A_, x_]:=(
    If[BuTools`CheckInput && Not[CheckMERepresentation[alpha, A]],Throw["CdfFromME: Input is not a valid ME distribution!"]];
	If[ListQ[x],
		Return[Table[1 - Total[alpha.MatrixExp[A xv]],{xv,x}]]
	, Return[1 - Total[alpha.MatrixExp[A x]]]]);


CdfFromPH[alpha_, A_, x_]:=(
    If[BuTools`CheckInput && Not[CheckPHRepresentation[alpha, A]],Throw["CdfFromPH: Input is not a valid PH distribution!"]];
	Return[CdfFromME[alpha,A,x]]);


CheckMEPositiveDensity[alpha_, A_, maxSize_:100, prec_:N[10^-12]]:=
Module[{beta,B},   
    Catch[
        {beta, B} = MonocyclicPHFromME[alpha, A, maxSize, prec];
        Return[CheckMERepresentation[beta, B, prec]]];
    Return[False]
];


CheckMERepresentation[alpha_, A_, precision_:Null]:=
Module[{size,eig,maxev,prec},

If[Not[NumericQ[precision]], prec=BuTools`CheckPrecision, prec=precision];

If[ Dimensions[A][[1]]!=Dimensions[A][[2]],
If[BuTools`Verbose,Print["CheckMERepresentation: The matrix is not a square matrix!"]];
Return[False]];

If[ Dimensions[alpha][[1]]!=Dimensions[A][[1]],
If[BuTools`Verbose,Print["CheckMERepresentation: The vector and the matrix have different sizes!"]];
Return[False]];

If[ Total[alpha]<-prec Length[alpha] || Total[alpha]>1+prec Length[alpha],
If[BuTools`Verbose,Print["CheckMERepresentation: The sum of the vector elements is not 1 (precision:",prec,")!"]];
Return[False]];

eig=Eigenvalues[A//N];
If[Max[Re[eig]]>=prec,
If[BuTools`Verbose,Print["CheckMERepresentation: There is an eigenvalue of the matrix with not negative real part at precision ",prec,")!"]];
 Return[False]];

eig=Sort[eig,Abs[Re[#1]]<Abs[Re[#2]] &];
maxev=eig[[1]];

If[Abs[Im[maxev]]>prec,
If[BuTools`Verbose,Print["CheckMERepresentation: The dominant eigenvalue of the matrix is not real at precision ",prec,")!"]];
 Return[False]];

If[Length[Select[eig,Abs[Im[#]]>prec && Re[#]==Re[eig[[1]]]&]]>0 && BuTools`Verbose,
Print["CheckMERepresentation: The dominant and a complex eigenvalue has the same real part at precision ",prec,")!!! "]];

Return[True] ;
];


CheckPHRepresentation[alpha_, A_,precision_:Null]:=
Module[{prec},
If[Not[NumericQ[precision]], prec=BuTools`CheckPrecision, prec=precision];
If[Dimensions[alpha][[1]]!=Dimensions[A][[1]],
If[BuTools`Verbose==True,Print["CheckPHRepresentation: the vector and the matrix have different sizes!"]];
Return[False]];
If[Not[CheckProbVector[alpha,True,prec]],Return[False],Return[CheckGenerator[A,True,prec]]];
];


IntervalPdfFromPH[alpha_, A_, intBounds_]:=
Module[{K,x},
    If[BuTools`CheckInput && Not[CheckPHRepresentation[alpha, A]],Throw["IntervalPdfFromPH: Input is not a valid PH distribution!"]];
    K = Length[intBounds];
	x = (intBounds[[2;;]] + intBounds[[1;;-2]])/2;
    Return[{x,Table[(Total[alpha.MatrixExp[A intBounds[[i]]]]-Total[alpha.MatrixExp[A intBounds[[i+1]]]])/(intBounds[[i+1]]-intBounds[[i]]),{i,1,K-1}]}];
];


MEFromMoments[moms_]:=
Module[{NN,KK,TT,UU,alpha,A,Appie},

	Appie[rmom_] :=
	Module[ {rm,m, i,j,f,y,dd,nold,yold,k,q,d,\[Alpha],\[Beta],n,\[Rho],kk,ind,inc},
		m = Length[rmom];
		If[Mod[m,2]==0, rm=Drop[rmom,-1];m=m/2, rm=rmom;m=Ceiling[m/2]];
		rm=Prepend[rm,1];
		f=ConstantArray[0,{2m,1}];
		f[[1,1]]=1;
		y=ConstantArray[0,{2m,1}];
		dd=ConstantArray[0,{2m,2m}];
		n=0;
		k=0;
		q=1;
		d=ConstantArray[0,{m}];
		\[Alpha]=ConstantArray[0,{m,m}];
		\[Beta]=ConstantArray[0,{m}];
		For[i=2,i<=2m,dd[[i,i-1]]=1;i++];
		For[i=1,i<=2m,
			\[Rho]=FullSimplify[q (rm.f)[[1]] ];
			nold=n;
			n=nold+1;
			yold=y;
			If[n>0 &&\[Rho]!=0,
				If[k>0,\[Beta][[k]]=\[Rho]/(rm[[1+1]])^(d[[k]]+n-1)];
				k=k+1;
				d[[k]]=n;
				n=-n;
				q=q/\[Rho];
				y=dd.f;
			, If[n<=0,
				j=nold+d[[k]]+1;
				\[Alpha][[k,j]]=\[Rho]/(rm[[1+1]])^(j-1);
			], (*UNDEFINED /in case of symbolic analysis/ *)
			If[n>0 ,
				If[k>0,\[Beta][[k]]=\[Rho]/(rm[[1+1]])^(d[[k]]+n-1)];
				k=k+1;
				d[[k]]=n;
				n=-n;
				q=q/\[Rho];
				y=dd.f;
			, If[n<=0,
				j=nold+d[[k]]+1;
				\[Alpha][[k,j]]=\[Rho]/(rm[[1+1]])^(j-1);
			]]];
			f= dd.f-\[Rho] yold;
			i++;
		];
		If[Sum[d[[i]],{i,1,m}]!=m,Print["Insufficient matrix order !!!!"]];
		kk=ConstantArray[0,{m,m}];
		kk[[1,1]]=rm[[2]];
		For[i=1,i<=m-1,kk[[i,i+1]]=rm[[2]];i++];
		ind=d[[1]];
		For[i=2,i<=m,
			If[ind<m,
				inc=d[[i]];
				If[(ind=ind+inc)<=m,
					kk[[ind,ind-inc-d[[i-1]]+1]]=\[Beta][[i-1]];
					For[j=1,j<=inc,
						kk[[ind,ind-j+1]]=\[Alpha][[i,j]];
						j++];
				];
			];
		i++];
		Return[kk];
	];

	KK = Appie[ReducedMomsFromMoms[moms]];
	NN = Ceiling[Length[moms]/2];

	TT = ConstantArray[0,{NN,NN}];
	Do[TT[[i,;;i]]=1,{i,NN}];

	UU = ConstantArray[0,{NN,NN}];
	Do[UU[[i,i;;]]=1/(NN-i+1),{i,NN}];

	alpha = ConstantArray[0,{NN}];
	alpha[[1]] = 1;
	alpha = alpha.Inverse[TT].UU;
	A = Inverse[-Inverse[UU].TT.KK.Inverse[TT].UU];
	Return[{alpha,A}];
] 


MEOrder[alpha_, A_, kind_:"moment", prec_:N[10^-10]]:=
Module[{NN,order,obsOrder,contOrder},
	If[BuTools`CheckInput && Not[CheckMERepresentation[alpha, A]],Throw["MinimalRepFromME: Input is not a valid ME distribution!"]];
	NN = Length[alpha];
	If[kind=="cont",
		order = MatrixRank[Table[Total[MatrixPower[Transpose[A],n-1],{1}],{n,NN}],Tolerance->prec];
    , If[kind=="obs",
		order = MatrixRank[Table[alpha.MatrixPower[A,n-1],{n,NN}],Tolerance->prec];
	, If[kind=="obscont",
		obsOrder = MatrixRank[Table[alpha.MatrixPower[A,n-1],{n,NN}],Tolerance->prec];
		contOrder = MatrixRank[Table[Total[MatrixPower[Transpose[A],n-1],{1}],{n,NN}],Tolerance->prec];
		order = Min[obsOrder,contOrder];
	, If[kind=="moment",
		order = MEOrderFromMoments[MomentsFromME[alpha, A], prec];
	, Throw["MEOrder: Invalid 'kind' parameter!"];
	]]]];
	Return[order];
];


MEOrderFromMoments[moms_,prec_ : N[10^-10] ]:=
Module[{size,rmoms,n,hankel,order},
	
	size = Floor[(Length[moms]+1)/2];
	rmoms = Prepend[ReducedMomsFromMoms[moms],1];
	order = size;
	For[n=1,n<=size,n++,
		hankel=Table[rmoms[[i+j-1]],{i,1,n},{j,1,n}];
		If[Abs[Det[hankel]] < prec, order=n-1; Break[]];
	];
	Return[order];
];


MinimalRepFromME[alpha_, A_, how_:"moment", precision_: N[10^-10]]:=
Module[{H0,H1,h,n,beta,B,alphav,Av,mo},
	If[BuTools`CheckInput && Not[CheckMERepresentation[alpha, A]],Throw["MinimalRepFromME: Input is not a valid ME distribution!"]];
	If[how=="cont",
		h = ConstantArray[1,{Length[alpha],1}];
		H0=A; 
		H1=-A.h.{alpha};
		{B,n} = MStaircase[{H0,H1},h,precision];
		beta = alpha.B;
		beta = beta[[1;;n]];
		B=Inverse[B].A.B;
		B=B[[1;;n,1;;n]];
    , If[how=="obs",
		h = ConstantArray[1,{Length[alpha],1}];
		H0=A; 
		H1=-A.h.{alpha};
		{B,n}=MStaircase[{Transpose[H0],Transpose[H1]},Transpose[{alpha}],precision];
		beta = alpha.B;
		beta = beta[[1;;n]];		
		B=Inverse[B].A.B;
		B=B[[1;;n,1;;n]];
	, If[how=="obscont",
		{alphav,Av}=MinimalRepFromME[alpha,A,"cont",precision];
		{beta,B}=MinimalRepFromME[alphav,Av,"obs",precision];
	, If[how=="moment",
        mo = MEOrder[alpha, A, "moment", precision];
        {beta, B} = MEFromMoments[MomentsFromME[alpha, A, 2 mo-1]];
	, Throw["MinimalRepFromME: Invalid 'how' parameter!"];
	]]]];
	Return[{beta,B}];
];


MomentsFromME[alpha_,A_,K_:0]:=
Module[{KK,iA},
    If[BuTools`CheckInput && Not[CheckMERepresentation[alpha, A]],Throw["MomentsFromME: Input is not a valid ME representation!"]];
    If[K==0,KK=2 Length[alpha]-1,KK=K];
	iA=Inverse[-A];
    Return[Table[i! Total[alpha.MatrixPower[iA,i]],{i,KK}]];
];


MomentsFromPH[alpha_,A_,K_:0]:=(
    If[BuTools`CheckInput && Not[CheckPHRepresentation[alpha, A]],Throw["MomentsFromPH: Input is not a valid PH representation!"]];
    MomentsFromME[alpha,A,K]);


MonocyclicPHFromME[alpha_, A_, maxSize_:100, precision_:N[10^-14]]:=
Module[{G,T,gamma,beta,B},
	If[BuTools`CheckInput && Not[CheckMERepresentation[alpha, A]],Throw["MonocyclicPHFromME: Input is not a valid ME distribution!"]];
    G = TransformToMonocyclic[A, maxSize, precision];
    T = SimilarityMatrix[A, G];
    gamma = Re[alpha.T];
    If[Min[gamma]>-precision,
		beta=gamma;B=G,
        {beta, B} = ExtendToMarkovian[gamma, G, maxSize, precision];
	];
    If[Not[CheckPHRepresentation[beta, B, precision]],Throw["MonocyclicPHFromME: No acyclic representation found up to the given size and precision!"]];
	Return[{beta,B}];
];


PdfFromME[alpha_, A_, x_]:=(
    If[BuTools`CheckInput && Not[CheckMERepresentation[alpha, A]],Throw["CdfFromME: Input is not a valid ME distribution!"]];
	If[ListQ[x],
		Return[Table[Total[alpha.MatrixExp[A xv].(-A)],{xv,x}]]
	, Return[Total[alpha.MatrixExp[A x].(-A)]]]);


PdfFromPH[alpha_, A_, x_]:=(
    If[BuTools`CheckInput && Not[CheckPHRepresentation[alpha, A]],Throw["PdfFromPH: Input is not a valid PH distribution!"]];
	Return[PdfFromME[alpha,A,x]]);


PH2From3Moments[moms_, prec_:N[10^-14]]:=
Module[{m1,m2,m3,m2l,m3l,m3u,a,b,c,e,lambda1,lambda2,p,alpha,A},
	{m1,m2,m3} = moms;
	(* check moment boounds*)
	m2l = APH2ndMomentLowerBound[m1, 2];
	m3l = APH3rdMomentLowerBound[m1, m2, 2];
	m3u = APH3rdMomentUpperBound[m1, m2, 2];
  
	If[m2<m2l, Throw["The given second moment is not feasible!"]];
	If[m3<m3l, Throw["The given third moment is not feasible (too small)!"]];
	If[m3>m3u, Throw["The given third moment is not feasible (too large)!"]];
    
	(* check if we have an exponential distribution*)
	If[Abs[m2/m1/m1-2]<prec, 
		alpha={1};
		A={{-1/m1}};
	,
		(* calculate parameters *)
		b = 3*m1*m2-m3;
		c = 3*m2*m2-2*m1*m3;
		e = -2*m1*m1+m2;
		a = b*b+6*c*e;
		If[a<0,a=0];
		a = Sqrt[a];
		If[c>0,
			lambda1 = (b - a) / c;
			lambda2 = (b + a) / c;
			p = (-b-6*m1*e+a) / (b+a);
		, If[c<0,
			lambda1 = (b + a) / c;
			lambda2 = (b - a) / c;
			p = (b+6*m1*e+a) / (-b+a);
		,
		    lambda1 = 0;
			lambda2 = 1 / m1;
			p = 0;
		]];
		alpha={p,1-p};
		A={{-lambda1,lambda1},{0,-lambda2}};
	];
	Return[{alpha,A}];
];


PH3From5Moments[moms_, prec_:N[10^-10]]:=
Module[{rmoms,M,m,a,discr,gu,g0,fs,ix,lambda,d1,d2,d3,gl,g2,
x1,x13,x2,x3,bels,p1,p2,p3,A,alpha},
    (* convert the moments to reduced moments*)
    rmoms = ReducedMomsFromMoms[moms];
	Do[rmoms[[i]]=rmoms[[i]]/moms[[1]]^i,{i,5}];
    (*solve linear system of equations for a0 a1 a2*)
    M = {{rmoms[[3]], -rmoms[[2]], rmoms[[1]]},{rmoms[[4]], -rmoms[[3]], rmoms[[2]]},{rmoms[[5]], -rmoms[[4]], rmoms[[3]]}};
	m = {1,rmoms[[1]],rmoms[[2]]};
	a = LinearSolve[M,m];    
    discr = a[[3]]^2-3*a[[2]];
    If[discr<0,Throw["Invalid characteristic polynomial!"]];
    
    gu = (a[[3]] + 2*Sqrt[discr]) / 3;
    g0 = (a[[3]] + Sqrt[discr]) / 3;

	fs=Map[#1[[1,2]]&,Solve[f^3+a[[3]] f^2+a[[2]] f+a[[1]]==0,f]];
    ix=Ordering[Re[fs]];
    lambda = -fs[[ix]];

    d1 = a[[2]] - a[[3]] - a[[1]] * rmoms[[2]];
    d2 = a[[1]] - a[[2]] - a[[3]] * d1;
    d3 = -a[[1]] - a[[2]]*d1 - a[[3]]*d2;

    If[d1>prec || (Abs[d1]<prec && d2>0), Throw["Negative density around 0!"]];
    If[lambda[[3]]<0,Throw["Invalid eigenvalues!"]];
    If[Im[lambda[[1]]] < prec, gl = Re[lambda[[1]]], gl = g0];
    If[gl>gu+prec, Throw["Invalid eigenvalues (gl>gu detected)!"]];
    If[gl>gu,gl = gu];
    If[Abs[d1]<prec, g2 = 0, g2 = -d2 / d1];
    If[g2>gu+prec,Throw["alpha_2 is negative!"]];
    If[g2>gu,g2 = gu];
    x1 = Max[g2, gl];
    If[Re[lambda[[1]]]==lambda[[1]] && g2<gl, x13 = 0, x13 = x1 - a[[1]] / (x1^2 - a[[3]]*x1 + a[[2]])];  
    bels = (a[[3]]-x1)^2 - 4*(x1^2-a[[3]]*x1+a[[2]]);
    If[bels<0 && bels>-prec, bels = 0];

    x2 = (a[[3]] - x1 + Sqrt[bels]) / 2;
    x3 = (a[[3]] - x1 - Sqrt[bels]) / 2;
    p1 = d1 / (x13 - x1);
    p2 = (x1*d1 + d2) / (x13-x1) / x2;
    p3 = (x1*x2*d1 + x2*d2 + x1*d2 + d3) / (x13-x1) / x2 / x3;

    A = {{-x1, 0, x13},{x2, -x2, 0},{0, x3, -x3}} / moms[[1]];
    alpha = {p1, p2, p3};

    If[x13<-prec || x13>x1, Throw["Invalid generator!"]];
    If[Min[Re[alpha]]<-prec, Throw["Initial vector has negative entries!"]];
    If[Max[Abs[Im[alpha]]]>prec,Throw["Inital vector has complex entries!"]];
    If[Max[Re[alpha]]>1+prec,Throw["Initial vector has entries that are greater than 1!"]];
    Return[{alpha,A}];
];


PHFromME[alpha_, A_, precision_:N[10^-14]]:=
Module[{Transfun,Evalfun,nrep},
    Transfun[orep_, B_]:=
        Return[{orep[[1]].B, Inverse[B].orep[[2]].B}];
  
    Evalfun[orep_, k_:0]:=Module[{ao,Ao,av,Ad},
        ao = orep[[1]];
        Ao = orep[[2]];
        av = Total[-Ao,{2}];
        Ad = Ao - DiagonalMatrix[Diagonal[Ao]];
        If[Mod[k,2] == 0,
            Return[-Min[ao, av, Ad]];
        ,
            Return[-Total[Select[ao,Negative]] - Total[Select[av,Negative]] - Total[Select[Flatten[Ad],Negative]]];
        ]];
   
    If[BuTools`CheckInput && Not[CheckMERepresentation[alpha, A]],Throw["PHFromME: Input is not a valid ME distribution!"]];
    nrep = FindMarkovianRepresentation[{alpha, A}, Transfun, Evalfun, precision];
    Return[nrep];
];


RandomPH[order_,mean_:1,zeroEntries_:0,maxTrials_:1000,prec_: N[10^-7]]:=
Module[{trials,zeros,numZeros,idx,zeroInRow,Row,B,sumRow,elem,vector,m,actualZeros},
If[zeroEntries>(order+1)(order-1),Throw["RandomPH: You have given too many zeros! Try to decrease the zero entries number!"]];
actualZeros=zeroEntries;
While[actualZeros>=0,
trials=1;
While[trials <= maxTrials,
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
If[zeros[[order+1]]>=order,trials++;Continue[]];(*Too many zeros in the last row*)
B=Table[0,{i,order},{j,order}];
For[idx=1,idx<=order,idx++,
Row=Table[Random[],{i,order-zeros[[idx]]}];
Row=Join[Row,Table[0,{i,zeros[[idx]]}]];
Row=Row[[RandomSample[Range[order],order]]];
B[[idx,1;;idx-1]]=Row[[1;;idx-1]];
B[[idx,idx+1;;order]]=Row[[idx;;order-1]];
B[[idx,idx]]=-Total[Row];
];
If[MatrixRank[B]!=order,trials++;Continue[]];
sumRow=0;
vector=Table[elem=Random[Real,1-sumRow];sumRow=sumRow+elem;elem,{i,order-zeros[[order+1]]-1}];
vector=Join[vector,{1-sumRow}];
vector=Join[vector,Table[0,{i,zeros[[order+1]]}]];
vector=vector[[RandomSample[Range[order],order]]];
If[Min[Abs[vector.Inverse[-B]]] < prec,trials++;Continue[]];

m=MomentsFromPH[vector,B,1];
B=B*m[[1]]/mean;
If[zeroEntries>actualZeros,
Print["RandomPH: Number of zero entries are different! Given:",zeroEntries,", in the returned representation:",actualZeros]];
Return[{vector,B}];
];
--actualZeros;
];
];


SamplesFromPH[a_,A_,k_]:=
Module[{NN,cummInitial,sojourn, nextpr,x,time,r,state,nstate},
	If[BuTools`CheckInput && Not[CheckPHRepresentation[a,A]],Throw["SamplesFromPH: input is not a valid PH representation!"]];
    (* auxilary variables*)
    NN = Length[a];
    cummInitial = Accumulate[a];
    sojourn = -1/Diagonal[A];
    nextpr = DiagonalMatrix[sojourn].A;
    nextpr = nextpr - DiagonalMatrix[Diagonal[nextpr]];
    nextpr = Join[nextpr, Transpose[{1-Total[nextpr,{2}]}],2];
    nextpr = Transpose[Accumulate[Transpose[nextpr]]];
    
    x = ConstantArray[0,{k}];
	Do[
		time = 0;
	    (* draw initial distribution *)
        r = RandomReal[];
        state = 1;
        While[cummInitial[[state]]<=r, state++];
        (* play state transitions *)
        While[state<=NN,
            time -= Log[RandomReal[]] sojourn[[state]];
            r = RandomReal[];
            nstate = 1;
            While[nextpr[[state,nstate]]<=r, nstate++];
            state = nstate;
        ];
        x[[n]] = time;
	,{n,k}];
	Return[x];
];


TestPHPackage[]:=
Module[{a,A,b,B,moms,err1,err2,x,cdf,pdf,y,flag,normMoms,phmoms,memoms,mean,n,mom2,cv2,
mp,Cm,co,oo,coo,mo,ma,mb,mom3lower,mom3upper},
Print["---BuTools: PH package test fil0^--"];
Print["Enable the verbose messages with the BuToolsVerbose flag"];
BuTools`Verbose = True;
Print["Enable input parameter checking with the BuToolsCheckInput flag"];
BuTools`CheckInput = True;
On[Assert];
Print["----------------------------"];
?MomentsFromME;

Print["Input:"];
Print["------"];

a = {0.2, 0.3, 0.5}; Print[a];
A = {{-1,0,0},{0,-3,2},{0,-2,-3}}; Print[A];

Print["Test:"];
Print["-----"];

Print["MomentsFromME[a,A]:"];
moms=MomentsFromME[a,A];
Print[moms];

Assert[CheckMoments[moms]==True, "The function returned invalid moments!"];

Print["Test:"];
Print["-----"];

Print["MomentsFromME[a,A,9]:"];
moms=MomentsFromME[a,A,9];
Print[moms];

Assert[CheckMoments[moms]==True, "The function returned invalid moments!"];

Print["----------------------------"];
?MomentsFromPH;

Print["Input:"];
Print["------"];

a={0.1, 0.9, 0}; Print[a];
A={{-6.2, 2, 0},{2, -9, 1},{1, 0, -3}}; Print[A];

Print["Test:"];
Print["-----"];

Print["MomentsFromPH[a,A,5)"];
moms=MomentsFromPH[a,A,5];
Print[moms];

Assert[CheckMoments[moms]==True, "The function returned invalid moments!"];

Print["----------------------------"];
?CdfFromME;

Print["Input:"];
Print["------"];

a = {0.2, 0.3, 0.5}; Print[a];
A = {{-1,0,0},{0,-3,2},{0,-2,-3}}; Print[A];

Print["Test:"];
Print["-----"];

Print["cdf=CdfFromME[a,A,Range[0,5,0.01]]:"];
cdf = CdfFromME[a,A,Range[0,5,0.01]];

Assert[AllTrue[Differences[cdf],Positive], "The cdf is not increasing monotonously!"];
Assert[Abs[Total[1-cdf]*0.01 - MomentsFromME[a,A,1][[1]]]<0.01, "The mean computed from the cdf does not match the theoretical result!"];

Print["----------------------------"];
?CdfFromPH;

Print["Input:"];
Print["------"];

a={0.1, 0.9, 0}; Print[a];
A={{-6.2, 2, 0},{2, -9, 1},{1, 0, -3}}; Print[A];

Print["Test:"];
Print["-----"];

Print["cdf=CdfFromPH[a,A,Range[0,3,0.002]]:"];
cdf = CdfFromPH[a,A,Range[0,3,0.002]];

Assert[AllTrue[Differences[cdf],Positive], "The cdf is not increasing monotonously!"];
Assert[Abs[Total[1-cdf]*0.002 - MomentsFromME[a,A,1][[1]]]<0.002, "The mean computed from the cdf does not match the theoretical result!"];

Print["----------------------------"];
?PdfFromME;

Print["Input:"];
Print["------"];

a = {0.2, 0.3, 0.5}; Print[a];
A = {{-1,0,0},{0,-3,2},{0,-2,-3}}; Print[A];
x = Range[0,5,0.01];

Print["Test:"];
Print["-----"];

Print["pdf=PdfFromME[a,A,x]:"];
pdf = PdfFromME[a,A,x];

Assert[AllTrue[pdf,NonNegative], "The pdf is negative!"];
Assert[Abs[pdf.x*0.01 - MomentsFromME[a,A,1][[1]]]<0.01, "The mean computed from the pdf does not match the theoretical result!"];

Print["----------------------------"];
?PdfFromPH;

Print["Input:"];
Print["------"];

a = {0.1, 0.9, 0}; Print[a];
A = {{-6.2, 2, 0},{2, -9, 1},{1, 0, -3}}; Print[A];
x = Range[0,3,0.002];

Print["Test:"];
Print["-----"];

Print["pdf=PdfFromPH[a,A,x]:"];
pdf = PdfFromPH[a,A,x];

Assert[AllTrue[pdf,NonNegative], "The pdf is negative!"];
Assert[Abs[pdf.x*0.002 - MomentsFromPH[a,A,1][[1]]]<0.002, "The mean computed from the pdf does not match the theoretical result!"];

Print["----------------------------"];
?IntervalPdfFromPH;

Print["Input:"];
Print["------"];

a = {0.1, 0.9, 0}; Print[a];
A = {{-6.2, 2, 0},{2, -9, 1},{1, 0, -3}}; Print[A];
x = Range[0,3,0.002];

Print["Test:"];
Print["-----"];

Print["{x,y}=IntervalPdfFromPH[a,A,x]:"];
{x,y} = IntervalPdfFromPH[a,A,x];

Assert[AllTrue[y,NonNegative], "The interval pdf is negative!"];
Assert[Abs[y.x*0.002 - MomentsFromPH[a,A,1][[1]]]<0.002, "The mean computed from the interval pdf does not match the theoretical result!"];

Print["----------------------------"];
?RandomPH;

Print["Test:"];
Print["-----"];
   
Print["[a,A]=RandomPH(3,8,4]:"];
{a,A}=RandomPH[3,8,4];
Print[a];
Print[A];

Assert[CheckPHRepresentation[a,A], "RandomPH failed to return a valid PH representation!"];
Assert[Max[Abs[MomentsFromPH[a,A,1][[1]]-8]]<10^-14, "RandomPH failed to match the given mean value!"];

Print["----------------------------"];
?CheckMERepresentation;

Print["Input:"];
Print["------"];

a={-0.2, 0.2}; Print[a];
A={{1, -1},{1, -2}}; Print[A];

Print["Test:"];
Print["-----"];

Print["CheckMERepresentation[a,A]:"];
flag=CheckMERepresentation[a,A];
Print[flag];

Assert[flag==False, "CheckMERepresentation did not detect that the initial vector is invalid!"];

Print["Input:"];
Print["------"];

a={-0.2, 0.4, 0.8}; Print[a];
A={{-2, 0, 3},{0, -1, 1},{ 0, -1, -1}}; Print[A];

Print["Test:"];
Print["-----"];

Print["CheckMERepresentation[a,A]:"];
flag=CheckMERepresentation[a,A];
Print[flag];

Assert[flag==False, "CheckMERepresentation did not detect that the dominant eigenvalue is invalid!"];

Print["Input:"];
Print["------"];

a = {0.2, 0.3, 0.5}; Print[a];
A = {{-1,0,0},{0,-3,2},{0,-2,-3}}; Print[A];

Print["Test:"];
Print["-----"];

Print["CheckMERepresentation[a,A]:"];
flag=CheckMERepresentation[a,A];
Print[flag];

Assert[flag==True, "CheckMERepresentation did not recognize that the given ME representation is valid!"];

Print["----------------------------"];
?CheckPHRepresentation;

Print["Input:"];
Print["------"];

a={0.2}; Print[a];
A={{-1, 1},{1, -2}}; Print[A];

Print["Test:"];
Print["-----"];

Print["CheckPHRepresentation[a,A]:"];
flag=CheckPHRepresentation[a,A];
Print[flag];

Assert[flag==False, "CheckPHRepresentation did not recognize the wrong input dimensions!"];

Print["Input:"];
Print["------"];

a={0.2, 0.7}; Print[a];
A={{-1, 1},{1, -2}}; Print[A];

Print["Test:"];
Print["-----"];

Print["CheckPHRepresentation[a,A]:"];
flag=CheckPHRepresentation[a,A];
Print[flag];

Assert[flag==True, "CheckPHRepresentation did not recognize that the given PH representation is valid!"];

Print["----------------------------"];
?CheckMEPositiveDensity;

Print["Input:"];
Print["------"];

a = {0.2, 0.3, 0.5}; Print[a];
A = {{-1,0,0},{0,-3,2},{0,-2,-3}}; Print[A];

Print["Test:"];
Print["-----"];

Print["CheckMEPositiveDensity[a,A]:"];
flag=CheckMEPositiveDensity[a,A//N];
Print[flag];

Assert[flag==True, "CheckMEPositiveDensity did not recognize that the given ME distribution has positive density!"];

Print["Input:"];
Print["------"];

a = {0.2, 0.3, 0.5}; Print[a];
A = {{-1,0,0},{0,-3,2.9},{0,-2.9,-3}}; Print[A];

Print["Test:"];
Print["-----"];

Print["CheckMEPositiveDensity[a,A]:"];
flag=CheckMEPositiveDensity[a,A];
Print[flag];

Assert[flag==False, "CheckMEPositiveDensity did not recognize that the given ME distribution does not have positive density!"];

Print["----------------------------"];
?APHFrom3Moments;

Print["Input:"];
Print["------"];

moms = {10,125,8400}; Print[moms];
normMoms = NormMomsFromMoms[moms]; Print[normMoms];

Print["Test:"];
Print["-----"];

Print["{a,A}=APHFrom3Moments[moms]:"];
Print["Length[a]"];
{a,A}=APHFrom3Moments[moms//N];
Print[a];
Print[A];
Print[Length[a]];

Print["MomentsFromPH[a,A,3]:"];
phmoms = MomentsFromPH[a,A,3];
Print[phmoms];

Assert[Norm[(phmoms-moms)/moms]<10^-12, "APHFrom3Moments failed to match the given moments!"];

Print["Input:"];
Print["------"];

moms = {10,525,31400}; Print[moms];
normMoms = NormMomsFromMoms[moms]; Print[normMoms];

Print["Test:"];
Print["-----"];

Print["{a,A}=APHFrom3Moments[moms]:"];
Print["Length[a]"];
{a,A}=APHFrom3Moments[moms//N];
Print[a];
Print[A];
Print[Length[a]];

Print["MomentsFromPH[a,A,3]:"];
phmoms = MomentsFromPH[a,A,3];
Print[phmoms];

Assert[Norm[(phmoms-moms)/moms]<10^-12, "APHFrom3Moments failed to match the given moments!"];


Print["----------------------------"];
?PH2From3Moments;

Print["Input:"];
Print["------"];

moms = {10,160,3500}; Print[moms];

Print["Test:"];
Print["-----"];

Print["{a,A}=PH2From3Moments[moms]"];
{a,A}=PH2From3Moments[moms//N];
Print[a];
Print[A];

Print["MomentsFromPH[a,A,3]:"];
phmoms = MomentsFromPH[a,A,3];
Print[phmoms];

Assert[Norm[(phmoms-moms)/moms]<10^-12, "PH2From3Moments failed to match the given moments!"];

Print["Input:"];
Print["------"];

moms = {10,260,13500}; Print[moms];

Print["Test:"];
Print["-----"];

Print["{a,A}=PH2From3Moments[moms]"];
{a,A}=PH2From3Moments[moms//N];
Print[a];
Print[A];

Print["MomentsFromPH[a,A,3]:"];
phmoms = MomentsFromPH[a,A,3];
Print[phmoms];

Assert[Norm[(phmoms-moms)/moms]<10^-12, "PH2From3Moments failed to match the given moments!"];

Print["----------------------------"];
?PH3From5Moments;

Print["Input:"];
Print["------"];

a={0.1, 0.9, 0}; Print[a];
A={{-6.2, 2, 0},{2, -9, 1},{1, 0, -3}}; Print[A];
moms = MomentsFromPH[a,A]; Print[moms];

Print["Test:"];
Print["-----"];

Print["{a,A}=PH3From5Moments[moms]"];
{a,A}=PH3From5Moments[moms]; 
Print[a];
Print[A];

Print["MomentsFromPH[a,A,5]:"];
phmoms = MomentsFromPH[a,A,5];
Print[phmoms];

Assert[Norm[(phmoms-moms)/moms]<10^-12, "PH3From5Moments failed to match the given moments!"];

Print["Input:"];
Print["------"];

a = {0.2, 0.3, 0.5}; Print[a];
A = {{-1,0,0},{0,-3,0.5},{0,-0.5,-3}}; Print[A];
moms = MomentsFromME[a,A]; Print[moms];

Print["Test:"];
Print["-----"];

Print["{a,A}=PH3From5Moments[moms]"];
{a,A}=PH3From5Moments[moms]; 
Print[a];
Print[A];

Print["MomentsFromPH[a,A,5]:"];
phmoms = MomentsFromPH[a,A,5];
Print[phmoms];

Assert[Norm[(phmoms-moms)/moms]<10^-12, "PH3From5Moments failed to match the given moments!"];

Print["----------------------------"];
?MEFromMoments;

Print["Input:"];
Print["------"];

a={0.1, 0.9, 0}; Print[a];
A={{-6.2, 2, 0},{2, -9, 1},{1, 0, -3}}; Print[A];
moms=MomentsFromPH[a,A,5]; Print[moms];

Print["Test:"];
Print["-----"];

Print["{a,A}=MEFromMoments[moms]:"];
{a,A}=MEFromMoments[moms];
Print[a];
Print[A];

Print["MomentsFromME[a,A,5]:"];
memoms = MomentsFromME[a,A,5];
Print[memoms];

Assert[Norm[(memoms-moms)/moms]<10^-12, "MEFromMoments failed to match the given moments!"];

Print["----------------------------"];
?APH2ndMomentLowerBound;

Print["Input:"];
Print["------"];

mean = 1.9;
n = 4;

Print["Test:"];
Print["-----"];

Print["mom2 = APH2ndMomentLowerBound[mean,n]:"];
mom2 = APH2ndMomentLowerBound[mean,n];
Print[mom2];

cv2 = mom2/mean^2-1;

Assert[Abs[cv2-1/n]<10^-14, "APH2ndMomentLowerBound did not give the expected result!"];


Print["----------------------------"];
?APH3rdMomentLowerBound;
?APH3rdMomentUpperBound;

Print["Input:"];
Print["------"];

mean = 1.9;
mom2 = 5;
n = 3;

Print["Test:"];
Print["-----"];

Print["mom3lower = APH3rdMomentLowerBound[mean,mom2,n]:"];
mom3lower = APH3rdMomentLowerBound[mean,mom2,n];
Print[mom3lower];
Print["mom3upper = APH3rdMomentUpperBound[mean,mom2,n]:"];
mom3upper = APH3rdMomentUpperBound[mean,mom2,n];
Print[mom3upper];

Assert[mom3upper>mom3lower, "Lower bound is larger than the upper bound!"];

Print["Input:"];
Print["------"];

mean = 1.9;
mom2 = 5;
n = 4;

Print["Test:"];
Print["-----"];

Print["mom3lower = APH3rdMomentLowerBound[mean,mom2,n]:"];
mom3lower = APH3rdMomentLowerBound[mean,mom2,n];
Print[mom3lower];
Print["mom3upper = APH3rdMomentUpperBound[mean,mom2,n]:"];
mom3upper = APH3rdMomentUpperBound[mean,mom2,n];
Print[mom3upper];

Assert[mom3upper>mom3lower, "Lower bound is larger than the upper bound!"];
Assert[mom3upper==Infinity, "Upper bound must be infinity with 4 phases!"];

Print["----------------------------"];
?CanonicalFromPH2;

Print["Input:"];
Print["------"];

a={0.12, 0.88}; Print[a];\[AliasDelimiter]
A={{-1.28, 0},{3.94, -3.94}}; Print[A];

Print["Test:"];
Print["-----"];

Print["{b,B}=CanonicalFromPH2[a,A]:"];
{b,B}=CanonicalFromPH2[a,A];
Print[b];
Print[B];

Cm=SimilarityMatrix[A,B];
err1 = Norm[A.Cm-Cm.B];
err2 = Norm[a.Cm-b];

Print["Transformation errors:"];
Print[err1];
Print[err2];

Assert[err1<10^-12 && err2<10^-12, "Transformation to canonical PH(2) failed!"];

Print["----------------------------"];
?CanonicalFromPH3;

Print["Input:"];
Print["------"];

a={0.1, 0.9, 0}; Print[a];
A={{-6.2, 2, 0},{2, -9, 1},{1, 0, -3}}; Print[A];

Print["Test:"];
Print["-----"];

Print["{b,B}=CanonicalFromPH3[a,A]:"];
{b,B}=CanonicalFromPH3[a,A];
Print[b];
Print[B];

Cm=SimilarityMatrix[A,B];
err1 = Norm[A.Cm-Cm.B];
err2 = Norm[a.Cm-b];

Print["Transformation errors:"];
Print[err1];
Print[err2];

Assert[err1<10^-12 && err2<10^-12, "Transformation to canonical PH(3) failed!"];

Print["----------------------------"];
?AcyclicPHFromME;

Print["Input:"];
Print["------"];

a={-0.4, 1.4, 0}; Print[a];
A={{-4, 1, 1},{0, -2, 1},{1, 0, -8}}; Print[A];

Print["Test:"];
Print["-----"];

Print["{b,B}=AcyclicPHFromME[a,A]:"];
{b,B}=AcyclicPHFromME[a,A//N];
Print[b];
Print[B];

Print["Moments of [a,A) and (b,B)"];
ma=MomentsFromME[a,A,5];
mb=MomentsFromME[b,B,5];
Print[ma];
Print[mb];

Assert[Norm[(ma-mb)/ma]<10^-7, "Transformation to acyclic representation failed!"];

Print["----------------------------"];
?MonocyclicPHFromME;

Print["Input:"];
Print["------"];

a = {0.2, 0.3, 0.5}; Print[a];
A = {{-1,0,0},{0,-3,2},{0,-2,-3}}; Print[A];

Print["Test:"];
Print["-----"];

Print["{b,B}=MonocyclicPHFromME[a,A]:"];
{b,B}=MonocyclicPHFromME[a,A//N];
Print[b];
Print[B];

Print["Moments of [a,A) and (b,B)"];
ma=MomentsFromME[a,A,5];
mb=MomentsFromME[b,B,5];
Print[ma];
Print[mb];

Assert[Norm[(ma-mb)/ma]<10^-7, "Transformation to monocyclic representation failed!"];

Print["----------------------------"];
?PHFromME;

Print["Input:"];
Print["------"];

a={-0.4, 1.4}; Print[a];
A={{-3.8, 2},{2, -9}}; Print[A];

Print["Test:"];
Print["-----"];

Print["CheckMERepresentation[a,A]:"];
flag=CheckMERepresentation[a,A];
Print[flag];

Print["CheckPHRepresentation[a,A]:"];
flag=CheckPHRepresentation[a,A];
Print[flag];

Print["{b,B}=PHFromME[a,A)"];
{b,B}=PHFromME[a,A];
Print[b];
Print[B];

Print["Check the obtained PH, CheckPHRepresentation(b,B]:"];
flag=CheckPHRepresentation[b,B];
Print[flag];

Cm=SimilarityMatrix[A,B];
err1 = Norm[A.Cm-Cm.B];
err2 = Norm[a.Cm-b];

Assert[flag && err1<10^-12 && err2<10^-12, "Transformation to PH failed!"];

Print["Input:"];
Print["------"];

a={-0.5, 1.5}; Print[a];
A={{-3.8, 2},{2, -9}}; Print[A];

Print["Test:"];
Print["-----"];

Print["CheckMERepresentation[a,A]:"];
flag=CheckMERepresentation[a,A];
Print[flag];

Print["CheckPHRepresentation[a,A]:"];
flag=CheckPHRepresentation[a,A];
Print[flag];

Print["{b,B}=PHFromME[a,A]"];
{b,B}=PHFromME[a,A];
Print[b];
Print[B];

Print["Check the obtained PH, CheckPHRepresentation(b,B]:"];
flag=CheckPHRepresentation[b,B];
Print[flag];

Cm=SimilarityMatrix[A,B];
err1 = Norm[A.Cm-Cm.B];
err2 = Norm[a.Cm-b];

Assert[flag && err1<10^-12 && err2<10^-12, "Transformation to PH failed!"];

Print["----------------------------"];
?MEOrder;

Print["Input:"];
Print["------"];

a={1, 1, 1, 1, 1, 1}/6; Print[a];
A={{-1, 0, 0, 0, 0, 0},{0.5, -2, 1, 0, 0, 0},{1, 0, -3, 1, 0, 0},{1, 0, 1, -4, 1, 0},{4, 0, 0, 0, -5, 0},{5, 0, 0, 0, 0, -6}}; Print[A];

Print["Test:"];
Print["-----"];

Print["co=MEOrder[a,A,\"cont\")"];
co=MEOrder[a,A//N,"cont"];
Print[co];
Print["oo=MEOrder[a,A,\"obs\")"];
oo=MEOrder[a,A//N,"obs"];
Print[oo];
Print["coo=MEOrder[a,A,\"obscont\")"];
coo=MEOrder[a,A//N,"obscont"];
Print[coo];
Print["coo=MEOrder[a,A,\"moment\")"];
mo=MEOrder[a,A//N,"moment"];
Print[mo];

Assert[co==2, "Wrong controllability order returned!"];
Assert[oo==6, "Wrong observability order returned!"];
Assert[coo==2, "The minimum of the controllability and observability order is wrong!"];
Assert[mo==2, "Wrong moment order returned!"];

Print["Input:"];
Print["------"];

a={2, 1}/3; Print[a];
A={{-1, 1},{0, -3}}; Print[A];

Print["Test:"];
Print["-----"];

Print["co=MEOrder[a,A,\"cont\")"];
co=MEOrder[a,A//N,"cont"];
Print[co];
Print["oo=MEOrder[a,A,\"obs\")"];
oo=MEOrder[a,A//N,"obs"];
Print[oo];
Print["coo=MEOrder[a,A,\"obscont\")"];
coo=MEOrder[a,A//N,"obscont"];
Print[coo];
Print["coo=MEOrder[a,A,\"moment\")"];
mo=MEOrder[a,A//N,"moment"];
Print[mo];

Assert[co==2, "Wrong controllability order returned!"];
Assert[oo==1, "Wrong observability order returned!"];
Assert[coo==1, "The minimum of the controllability and observability order is wrong!"];
Assert[mo==1, "Wrong moment order returned!"];

Print["Input:"];
Print["------"];

b = {0.2, 0.3, 0.5}; Print[b];
B = {{-1,0,0},{0,-3,1},{0,-1,-3}}; Print[B];
{a,A} = MonocyclicPHFromME[b,B//N];
Print[a];
Print[A];

Print["Test:"];
Print["-----"];

Print["co=MEOrder[a,A,\"cont\")"];
co=MEOrder[a,A//N,"cont"];
Print[co];
Print["oo=MEOrder[a,A,\"obs\")"];
oo=MEOrder[a,A//N,"obs"];
Print[oo];
Print["coo=MEOrder[a,A,\"obscont\")"];
coo=MEOrder[a,A//N,"obscont"];
Print[coo];
Print["coo=MEOrder[a,A,\"moment\")"];
mo=MEOrder[a,A//N,"moment"];
Print[mo];

Assert[co==9, "Wrong controllability order returned!"];
Assert[oo==3, "Wrong observability order returned!"];
Assert[coo==3, "The minimum of the controllability and observability order is wrong!"];
Assert[mo==3, "Wrong moment order returned!"];

Print["----------------------------"];
?MEOrderFromMoments;

Print["Input:"];
Print["------"];

a={0.1, 0.9, 0}; Print[a];
A={{-6.2, 2, 0},{2, -9, 1},{1, 0, -3}}; Print[A];

Print["Test:"];
Print["-----"];

Print["moms=MomentsFromME[a,A)"];
Print["mo=MEOrderFromMoments[moms]:"];
moms=MomentsFromME[a,A];
mo = MEOrderFromMoments[moms];
Print[mo];

Assert[mo==3, "Wrong moment order returned!"];

Print["Input:"];
Print["------"];

b = {0.2, 0.3, 0.5}; Print[b];
B = {{-1,0,0},{0,-3,2},{0,-2,-3}}; Print[B];
{a,A} = MonocyclicPHFromME[b,B//N];
Print[a];
Print[A];

Print["Test:"];
Print["-----"];

Print["moms=MomentsFromME[a,A]"];
Print["mo=MEOrderFromMoments[moms]:"];
moms=MomentsFromME[a,A];
mo = MEOrderFromMoments[moms];
Print[mo];

Assert[mo==3, "Wrong moment order returned!"];

Print["----------------------------"];
?MinimalRepFromME;

Print["Input:"];
Print["------"];

a={1, 1, 1, 1, 1, 1}/6; Print[a];
A={{-1, 0, 0, 0, 0, 0},{0.5, -2, 1, 0, 0, 0},{1, 0, -3, 1, 0, 0},{1, 0, 1, -4, 1, 0},{4, 0, 0, 0, -5, 0},{5, 0, 0, 0, 0, -6}}; Print[A];

Print["Test:"];
Print["-----"];

Print["{b,B}=MinimalRepFromME[a,A,\"cont\"]"];
{b,B}=MinimalRepFromME[a,A//N,"cont"];
Print[b];
Print[B];

Assert[Length[b]==2, "Non-minimal representation returned based on controllability!"];

Print["{b,B}=MinimalRepFromME[a,A,\"obs\")"];
{b,B}=MinimalRepFromME[a,A//N,"obs"];
Print[b];
Print[B];

Assert[Length[b]==6, "Non-minimal representation returned based on observability!"];

Print["{b,B}=MinimalRepFromME[a,A, \"obscont\"]"];
{b,B}=MinimalRepFromME[a,A//N,"obscont"];
Print[b];
Print[B];

Assert[Length[b]==2, "Non-minimal representation returned based on observability and controllability!"];

Print["{b,B}=MinimalRepFromME[a,A,\"moment\"]"];
{b,B}=MinimalRepFromME[a,A//N,"moment"];
Print[b];
Print[B];

Assert[Length[b]==2, "Non-minimal representation returned based on the moments!"];

Print["Input:"];
Print["------"];

a={2, 1}/3; Print[a];
A={{-1, 1},{0, -3}}; Print[A];

Print["Test:"];
Print["-----"];

Print["{b,B}=MinimalRepFromME[a,A,\"cont\"]"];
{b,B}=MinimalRepFromME[a,A//N,"cont"];
Print[b];
Print[B];

Assert[Length[b]==2, "Non-minimal representation returned based on controllability!"];

Print["{b,B}=MinimalRepFromME[a,A,\"obs\"]"];
{b,B}=MinimalRepFromME[a,A//N,"obs"];
Print[b];
Print[B];

Assert[Length[b]==1, "Non-minimal representation returned based on observability!"];

Print["{b,B}=MinimalRepFromME[a,A,\"obscont\"]"];
{b,B}=MinimalRepFromME[a,A//N,"obscont"];
Print[b];
Print[B];

Assert[Length[b]==1, "Non-minimal representation returned based on observability and controllability!"];

Print["{b,B}=MinimalRepFromME[a,A,\"moment\"]"];
{b,B}=MinimalRepFromME[a,A//N,"moment"];
Print[b];
Print[B];

Assert[Length[b]==1, "Non-minimal representation returned based on the moments!"];

Print["Input:"];
Print["------"];

b = {0.2, 0.3, 0.5}; Print[b];
B = {{-1,0,0},{0,-3,1},{0,-1,-3}}; Print[B];
{a,A} = MonocyclicPHFromME[b,B//N];
Print[a];
Print[A];

Print["{b,B}=MinimalRepFromME[a,A,\"cont\"]"];
{b,B}=MinimalRepFromME[a,A//N,"cont"];
Print[b];
Print[B];

Assert[Length[b]==Length[a], "Non-minimal representation returned based on controllability!"];

Print["{b,B}=MinimalRepFromME[a,A,\"obs\"]"];
{b,B}=MinimalRepFromME[a,A//N,"obs"];
Print[b];
Print[B];
(* check similarity*)
Cm=SimilarityMatrix[B,A];
err1 = Norm[B.Cm-Cm.A];
err2 = Norm[b.Cm-a];

Assert[Length[b]==3 && err1+err2<10^-12, "Non-minimal representation returned based on observability!"];

Print["{b,B}=MinimalRepFromME[a,A,\"obscont\")"];
{b,B}=MinimalRepFromME[a,A//N,"obscont"];
Print[b];
Print[B];
(* check similarity:*)
Cm=SimilarityMatrix[B,A];
err1 = Norm[B.Cm-Cm.A];
err2 = Norm[b.Cm-a];

Assert[Length[b]==3 && err1+err2<10^-12, "Non-minimal representation returned based on observability and controllability!"];

Print["{b,B}=MinimalRepFromME[a,A,\"moment\")"];
{b,B}=MinimalRepFromME[a,A//N,"moment"];
Print[b];
Print[B];
(* check similarity:*)
Cm=SimilarityMatrix[B,A];
err1 = Norm[B.Cm-Cm.A];
err2 = Norm[b.Cm-a];

Assert[Length[b]==3 && err1+err2<10^-12, "Non-minimal representation returned based on the moments!"];

Print["----------------------------"];
?SamplesFromPH;

Print["Input:"];
Print["------"];

a={0.1, 0.9, 0}; Print[a];
A={{-6.2, 2, 0},{2, -9, 1},{1, 0, -3}}; Print[A];

Print["Test:"];
Print["-----"];

Print["x=SamplesFromPH[a,A,1000]"];
x=SamplesFromPH[a,A,1000];

(*Print["Moments from the samples:"];
mt = MarginalMomentsFromTrace[x,3];
Print[mt];
*)
Print["Moments from the PH:"];
mp = MomentsFromPH[a,A,3];
Print[mp];
];


End[(* Private *)];


EndPackage[];
