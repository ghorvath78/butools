(* ::Package:: *)

(*
   BuTools_PH Package
*)

(*FONTOS! A RUN PACKAGE CSAK AKKOR FOG JOL MUKODNI, HA ELOTTE*)
(*A PROBA FAJLBAN A DIR, APPENTO, PATH RESZT LEFUTTATOD!!! *)

BeginPackage["BuTools`PH`"];
RandomPH::usage = "RandomPH [ order, zeroEntries, mean[1], \[Epsilon][10^-14], maxTrials[100] ] -> [ vector, matrix ] : Generates a random continous time Phase Type process of the given order. The obtained representation containes 'zeroEntries' zeros, its first moment is 'mean'. If it fails after 'maxTrials' trials, then it decreases the number of zero entries. It prints a message, if the found representation contains less zeros than given. \[Epsilon] is the numerical precision.";
RandomDPH::usage = "RandomDPH [ order, zeroEntries, \[Epsilon][10^-14], maxTrials[100] ] -> [ vector, matrix ] : Generates a random continous time Phase Type process of the given order. The obtained representation containes 'zeroEntries' zeros. If it fails after 'maxTrials' trials, then it decreases the number of zero entries. It prints a message, if the found representation contains less zeros than given. \[Epsilon] is the numerical precision.";
PHFromME::usage = "PHFromME [ vector, matrix ] -> [ vector, matrix ] : Converts  matrix exponential vector-matrix pair to phase type vector-matrix pair.";
MomentsFromME::usage = "MomentsFromME [ vector, matrix, k[2n-1], \[Epsilon][10^-14] ] -> [ moments ] : Calculates the first k moments of a ME given with a vector-matrix pair. \[Epsilon] is the numerical precision.";
MomentsFromMG::usage = "MomentsFromMG [ vector, matrix, k[2n-1], \[Epsilon][10^-14] ] -> [ moments ] : Calculates the first k moments of a MG given with a vector-matrix pair. \[Epsilon] is the numerical precision.";
MomentsFromPH::usage = "MomentsFromPH [ vector, matrix, k[2n-1], \[Epsilon][10^-14] ] -> [ moments ] : Calculates the first k moments of a PH given with a vector-matrix pair. \[Epsilon] is the numerical precision.";
MomentsFromDPH::usage = "MomentsFromDPH [ vector, matrix, k[2n-1], \[Epsilon][10^-14] ] -> [ moments ] : Calculates the first k moments of a DPH given with a vector-matrix pair. \[Epsilon] is the numerical precision.";
MEFromMoments::usage = "MEFromMoments [ moments ] -> [ vector, matrix ] : Calculates the minimal ME representation based on 2n-1 moments.";
MGFromMoments::usage = "MGFromMoments [ moments ] -> [ vector, matrix ] : Calculates the minimal MG representation based on 2n-1 moments.";
PH2Canonical::usage = "PH2Canonical [ vector, matrix, \[Epsilon][10^-14] ] -> [ vector, matrix ] : Calculates the order 2 canonical representation from any order 2 vector-matrix representation, if exists. \[Epsilon] is the numerical precision.";
DPH2Canonical::usage = "DPH2Canonical [ vector, matrix, \[Epsilon][10^-14] ] -> [ vector, matrix ] : Calculates the order 2 canonical representation from any order 2 vector-matrix representation, if exists. \[Epsilon] is the numerical precision.";
PH3Canonical::usage = "PH3Canonical [ vector, matrix, \[Epsilon][10^-14] ] -> [ vector, matrix ] : Calculates the order 3 canonical representation from any order 3 PH representation. \[Epsilon] is the numerical precision.";
DPH3Canonical::usage = "DPH3Canonical [ vector, matrix, \[Epsilon][10^-14] ] -> [ vector, matrix ] : Calculates the order 3 canonical representation from any order 3 DPH representation. \[Epsilon] is the numerical precision.";
APHRepresentation::usage = "APHRepresentation [ vector, matrix, \[Epsilon][10^-14] ] -> [ vector, matrix ] : Calculates the APH (CF1) representation from any order n vector-matrix representation, if exists. \[Epsilon] is the numerical precision.";
ADPHRepresentation::usage = "ADPHRepresentation [ vector, matrix, \[Epsilon][10^-14] ] -> [ vector, matrix ] : Calculates the ADPH (CF1) representation from any order n vector-matrix representation, if exists. \[Epsilon] is the numerical precision.";
MonocyclicRepresentation::usage = "MonocyclicRepresentation [ vector, matrix, exitProb[0] ] -> [ vector, matrix ] : Calculates the representation of the input ME distribution with Markovian monocyclic generator.  The process terminates with 'exitProb' probability and jumps to the next FE block with 1-exitProb probability at every arrival.";
RepTrafo::usage = "RepTrafo [ inVector, inMatrix, outMatrix ] -> [ vector ] : Calculates the transformation from inMatrix to outMatrix, and returns the transformed inVector.";
Appie::usage = "Appie [ reducedmoments, print flag ] -> [ matrix ] : Calculates the minimal ME representation with e_1 opening and closing vector based on 2n-1 reducedmoments.";
APHFrom3Moments::usage = "APHFrom3Moments [ mom1, mom2, mom3] -> [ vector, matrix ] : Calculates the smallest APH with these 3 moments.";
MEOrderFromMoments::usage = "MEOrderFromMoments [ moments ] -> [ order ] : Calculates the order of the ME distribution based on its moments using the determinant of the Hankel matrix.";
ME3member::usage = "ME3member [ vector, matrix ] -> [ flag ] : Checks if the vector-matrix pair is a member of ME3.";
MEContOrder::usage = "MEContOrder [ vector, matrix ] -> [ order ] : Controllability (closing vector) order of the vector-matrix pair.";
MEObsOrder::usage = "MEObsOrder [ vector, matrix ] -> [ order ] : Observability (initial vector) order of the vector-matrix pair.";
MEContMinimize::usage = "MEContMinimize [ vector, matrix ] -> [ vector, matrix ] : Minimizies the representation according to the closing vector.";
MEObsMinimize::usage = "MEObsminimize [ vector, matrix ] -> [ vector, matrix ] : Minimizies the representation according to the initial vector.";
MEMinimize::usage = "MEMiminize [ vector, matrix ] -> [ vector, matrix ] : Minimizies the representation according to the initial and the closing vector.";
CheckMEPositiveDensity::usage = "CheckMEPositiveDensity [ vector, matrix ] -> [ flag ] : Checks if the vector-matrix pair results in a positive density.";
MEDensity::usage = "MEDensity [ vector, matrix, x ]  -> [ density value ] : Gives back the value of the probability density function at point x.";


Begin["`Private`"];


Needs["BuTools`RepTrans"];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckInput"]],BuTools`CheckInput=True];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckPrecision"]],BuTools`CheckPrecision=N[10^-12]];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`Verbose"]],BuTools`Verbose=False];


AcyclicPHFromME[alpha_, A_, maxSize_, precision_]:=
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
Module[{cv2,lambda,p,n},
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
Module[{size,eig,ceig,reig,prec},

If[Not[NumericQ[precision]], prec=BuTools`CheckPrecision, prec=precision];

If[ Dimensions[A][[1]]!=Dimensions[A][[2]],
If[BuTools`Verbose,Print["CheckMERepresentation: The matrix is not a square matrix!"]];
Return[False]];

If[ Dimensions[alpha][[1]]!=Dimensions[A][[1]],
If[BuTools`Verbose,Print["CheckMERepresentation: The vector and the matrix have different sizes!"]];
Return[False]];

If[ Total[alpha]<-prec Length[alpha] || Total[alpha]>1+prec Length[alpha],
If[BuTools`Verbose,Print["CheckMERepresentation: The sum of the vector elements is not 1 (precision:",\[Epsilon],")!"]];
Return[False]];

eig=Eigenvalues[A//N];
If[Max[Re[eig]]>=prec,
If[BuTools`Verbose,Print["CheckMERepresentation: There is an eigenvalue of the matrix with not negative real part at precision ",\[Epsilon],")!"]];
 Return[False]];

size=Dimensions[A][[1]];
ceig=Table[Switch[Element[eig[[i]],Reals],True,-Infinity,False,Re[eig[[i]]]],{i,size}];
reig=Table[Switch[Element[eig[[i]],Reals],True,Re[eig[[i]]],False,-Infinity],{i,size}];

If[Max[reig]+prec<Max[ceig],
If[BuTools`Verbose,Print["CheckMERepresentation: The dominant eigenvalue of the matrix is not real at precision ",\[Epsilon],")!"]];
 Return[False]];

If[Abs[Max[reig]-Max[ceig]]<prec && BuTools`Verbose,
Print["CheckMERepresentation: The dominant and a complex eigenvalue has the same real part at precision ",\[Epsilon],")!!! "]];

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


RandomPH[order_,zeroEntries_,mean_: 1,\[Epsilon]_: N[10^-14],maxTrials_: 100]:=
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
If[Min[Abs[vector.Inverse[-B]]] < \[Epsilon],trials++;Continue[]];

m=MomentsFromPH[vector,B,1,\[Epsilon]];
B=B*m[[1]]/mean;
If[zeroEntries>actualZeros,
Print["RandomPH: Number of zero entries are different! Given:",zeroEntries,", in the returned representation:",actualZeros]];
Return[{vector,B}];
];
--actualZeros;
];
];


MEFromMoments[mom_]:=
Module[{size,kk,rmom,hh,e1,tt,uu,i,j,tti,uui},
rmom=ReducedmomsFromMoms[mom];
{kk,size}=Appie[rmom,False] ;

hh=ConstantArray[1,{size,1}];
e1=ConstantArray[0,{1,size}];
e1[[1,1]]=1;

tt=ConstantArray[0,{size,size}];
For[i=1,i<=size,i++,For[j=1,j<=i,j++,tt[[i,j]]= 1;];];
tti=Inverse[tt];

uu=ConstantArray[0,{size,size}];
For[i=1,i<=size,i++,For[j=i,j<=size,j++,uu[[i,j]]= 1/(size-i+1);];];
uui=Inverse[uu];

Return[{(e1)[[1]].tti.uu,-Inverse[uui.tt.kk.tti.uu]}];
] 


PH3Canonical[v_,h0_,\[Epsilon]_ : N[10^-14] ] :=
Module[ {ii,\[Lambda],a0,a1,a2,a3,x1,x2,x3,x13,\[Gamma]u,\[Gamma]0,\[Gamma]l,\[Gamma]2,f,p1,p2,p3,p,e,\[Alpha]},

If[Not[CheckMERepresentation[v, h0, \[Epsilon]]], Throw["PH3Canonical: Input isn't a valid ME distribution!"]];
If[Dimensions[h0][[1]] !=3, Throw["PH3Canonical: Dimension \[NotEqual] 3!"] ];

If[Not[CheckPHRepresentation[v,h0,\[Epsilon]]], 
If[BuToolsVerbose,Print["PH3Canonical Warning: The input isn't a PH representation!"]];
];


e={{1},{1},{1}};

If[-(v.h0.e)[[1]]<-\[Epsilon],Throw["PH3Canonical: Negative density at 0!"]];
ii=IdentityMatrix[3];
\[Lambda]=Sort[Eigenvalues[-h0],Abs[#1]>Abs[#2]||(Abs[#1]==Abs[#2] && Re[#1]>Re[#2])||(Abs[#1]==Abs[#2]&&Re[#1]==Re[#2]&&Im[#1]>=Im[#2])&];
If[BuToolsVerbose,Print["eigenvalues: ",\[Lambda]]];

If[(Re[\[Lambda][[3]]]<0 )||Abs[Im[\[Lambda][[3]]]]>\[Epsilon], 
 Throw["PH3Canonical: eigenvalue error!!! \[Lambda]3=",\[Lambda][[3]]];
];

a0=\[Lambda][[1]] \[Lambda][[2]] \[Lambda][[3]] ;
a1=\[Lambda][[1]] \[Lambda][[2]] + \[Lambda][[1]] \[Lambda][[3]] + \[Lambda][[2]] \[Lambda][[3]];
a2=\[Lambda][[1]] + \[Lambda][[2]] +  \[Lambda][[3]];
{a0,a1,a2}=Re[{a0,a1,a2}];

f[x_]:= x^3+a2  x^2+a1 x+a0;
\[Gamma]u=(a2+2 Sqrt[a2^2-3 a1])/3;
\[Gamma]0=(a2+ Sqrt[a2^2-3 a1])/3;
If[\[Lambda][[1]]\[Element] Reals,\[Gamma]l=\[Lambda][[1]],\[Gamma]l=\[Gamma]0];

If[Abs[(v.h0.e)[[1]]]<\[Epsilon],\[Gamma]2=0,\[Gamma]2= -(v.h0.h0.e)[[1]]/(v.h0.e)[[1]]];

If[BuToolsVerbose,Print["\[Gamma]u=",\[Gamma]u,", \[Gamma]l=",\[Gamma]l,", \[Gamma]2=",\[Gamma]2]];
If[Max[\[Gamma]l,\[Gamma]2] > \[Gamma]u, If[BuToolsVerbose,Print["Max(\[Gamma]l,\[Gamma]2) > \[Gamma]u !!!!!!"]]];
x1=Max[\[Gamma]l,\[Gamma]2];
If[BuToolsVerbose,Print["Degree of freedom: \[Gamma]u - Max(\[Gamma]l,lower)=",\[Gamma]u - Max[\[Gamma]l,\[Gamma]2]]]; 

If[x1==\[Lambda][[1]],x13=0,x13=-f[-x1]/(a0-f[-x1]) x1];
x2=(a2-x1+Sqrt[(a2-x1)^2-4(x1^2 - a2 x1 +a1)])/2;
x3=(a2-x1-Sqrt[(a2-x1)^2-4(x1^2 - a2 x1 +a1)])/2;
If[BuToolsVerbose,Print["x1=",x1,", x13=",x13,", x2=",x2,", x3=",x3]];
{x13,x2,x3}=Re[{x13,x2,x3}];
If[x13>x1+\[Epsilon],Print["x13 > x1 !!!! x1=",x1,", x13=",x13] ];
If[x13<-\[Epsilon],Print["x13 < 0 !!!! x13=",x13] ];
p1=- h0.{{1},{1},{1}}/(x1-x13);p1=Transpose[p1][[1]];
p2=-(x1 ii + h0).h0.{{1},{1},{1}}/(x1-x13)/x2;p2=Transpose[p2][[1]];
p3=-(x2 ii + h0).(x1 ii + h0).h0.{{1},{1},{1}}/(x1-x13)/x2/x3;p3=Transpose[p3][[1]];
p=Transpose[{p1,p2,p3}];
\[Alpha]=v.p;
If[x1==\[Gamma]2,\[Alpha][[2]]=0];
Return[{\[Alpha],{{-x1,0,x13},{x2,-x2,0},{0,x3,-x3}}}];
];


Appie[rmom_,print_] :=
Module[ {rm,m, i,j,f,y,dd,nold,yold,k,q,d,\[Alpha],\[Beta],n,\[Rho],kk,ind,inc},
     m=Dimensions[rmom][[1]];
If[Mod[m,2]==0,
Print["Even number of moments, the last one is dropped !!"];
rm=Drop[rmom,-1];m=m/2,rm=rmom;m=Ceiling[m/2]
];
rm=Prepend[rm,1];
f=Table[0,{2m},{1}];
f[[1,1]]=1;
y=Table[0,{2m},{1}];
dd=Table[0,{2m},{2m}];
n=0;
k=0;
q=1;
d=Table[0,{m}];
\[Alpha]=Table[0,{m},{m}];
\[Beta]=Table[0,{m}];
For[i=2,i<=2m,dd[[i,i-1]]=1;i++];
For[i=1,i<=2m,
(*THEN*)
\[Rho]=FullSimplify[q (rm.f)[[1]] ];
nold=n;
n=nold+1;
yold=y;
If[print,Print["iter ",i," \[Rho]=",\[Rho]," n=",n, " yold=",yold]];
If[n>0 &&\[Rho]!=0 ,
(*THEN*)
If[k>0,\[Beta][[k]]=\[Rho]/(rm[[1+1]])^(d[[k]]+n-1)];
k=k+1;
d[[k]]=n;
n=-n;
q=q/\[Rho];
y=dd.f;
If[print,Print[i," then: \[Beta]=",\[Beta]," k=",k," n=",n," q=",q ]],
(*ELSE*)
If[print,Print[i," else: n=",n]];
If[n<=0,
j=nold+d[[k]]+1;
\[Alpha][[k,j]]=\[Rho]/(rm[[1+1]])^(j-1);
If[print,Print[i ," else if: k+1=",k+1," j+1=",j+1, " \[Alpha]=",\[Alpha]]]
],
(*UNDEFINED /in case of symbolic analysis/ *)
If[n>0 ,
(*THEN*)
If[k>0,\[Beta][[k]]=\[Rho]/(rm[[1+1]])^(d[[k]]+n-1)];
k=k+1;
d[[k]]=n;
n=-n;
q=q/\[Rho];
y=dd.f;
If[print,Print[i," undefined then: \[Beta]=",\[Beta]," k=",k," n=",n," q=",q ]],
(*ELSE*)
If[print,Print[i," else: n=",n]];
If[n<=0,
j=nold+d[[k]]+1;
\[Alpha][[k,j]]=\[Rho]/(rm[[1+1]])^(j-1);
If[print,Print[i ," undefined else if: k+1=",k+1," j+1=",j+1, " \[Alpha]=",\[Alpha]]]
]
]
];
f= dd.f-\[Rho] yold;
If[print,Print["f=",f]];
i++];

If[print,Print["m=",m," \[Beta]=",\[Beta]," \[Alpha]=",\[Alpha]," d=",d]];
If[Sum[d[[i]],{i,1,m}]!=m,Print["Insufficient matrix order !!!!"]];

kk=Table[0,{m},{m}];
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
Return[{kk,m}];
]


MEOrderFromMoments[moms_,\[Epsilon]_ : N[10^-14] ]:=
Module[{momsize,size,rmoms,i,j,n,m,hankel},
momsize= Dimensions[moms][[1]];
size= Floor[(momsize+1)/2];
rmoms=ReducedmomsFromMoms[moms];
rmoms=Prepend[rmoms,1];
For[n=1,n<=size,n++,
hankel=Table[rmoms[[i+j-1]],{i,1,n},{j,1,n}];
Print["n=",n,", det hankel=",Det[hankel], ", precision=",\[Epsilon]//N];
If[ (* Det[hankel]==0 *) Abs[Det[hankel]] < \[Epsilon] && n==size, Return[n-1]];

If[(* Det[hankel]==0 *) Abs[Det[hankel]] <\[Epsilon],
For[m=1,m<= momsize- 2n+1,m++,
hankel=Table[rmoms[[i+j-1+m]],{i,1,n},{j,1,n}];
If[(* Det[hankel]!=0 *) Abs[Det[hankel]] > \[Epsilon],Print["Hankel matrix of r(0)-r(",2n-1,") is 0, but the Hankel matrix of r(",m,")-r(",2n-1+m,") is not 0 !!"]; Return[n];
];
]; (* end For *)
Return[n-1];
];
];
Return[size];
];


MEContOrder[\[Alpha]_,A_]:=
Module[{size,h,n,res},
size=Dimensions[A][[1]];
h=Table[1,{size}];
res=MatrixRank[Table[h.MatrixPower[Transpose[A],n-1],{n,size}]];
Return[res];
];


MEObsOrder[\[Alpha]_,A_]:=
Module[{size,n,res},
size=Dimensions[A][[1]];
res=MatrixRank[Table[\[Alpha].MatrixPower[A,n-1],{n,size}]];
Return[res];
];


MEContMinimize[\[Alpha]in_,Ain_]:= 
Module[ {d0,d1,size,vec,B,n,i,rowsum},
size=Dimensions[Ain][[1]];
vec= Table[1,{size},{1}];
d0=Ain; d1=-Ain.vec.{\[Alpha]in};

{n,B}=StairCase[d0,d1,vec];
d0=Inverse[B].d0.B;
d1=Inverse[B].d1.B;
d0=d0[[1;;n,1;;n]];
d1=d1[[1;;n,1;;n]];

For[i=1,i<=n,i++,
rowsum=d1[[i]].Table[1,{n}];
If[rowsum!=0,
Return[{d1[[i]]/rowsum,d0}]]
];
Return[{d1[[1]],d0}];
]


MEObsMinimize[\[Alpha]in_,Ain_]:= 
Module[ {d0,d1,size,vec,B,n,i,rowsum},
size=Dimensions[Ain][[1]];
vec= Table[1,{size},{1}];
d0=Ain; d1=-Ain.vec.{\[Alpha]in};

vec={\[Alpha]in};
{n,B}=StairCase[Transpose[d0],Transpose[d1],Transpose[vec]];
d0=Inverse[B].d0.B;
d1=Inverse[B].d1.B;
d0=d0[[1;;n,1;;n]];
d1=d1[[1;;n,1;;n]];

For[i=1,i<=n,i++,
rowsum=d1[[i]].Table[1,{n}];
If[rowsum!=0,
Return[{d1[[i]]/rowsum,d0}]]
];
Return[{d1[[1]],d0}];

]


MEMinimize[\[Alpha]in_,Ain_]:= 
Module[ {\[Alpha], A},
{\[Alpha], A}=MEContMinimize[\[Alpha]in,Ain]//Chop ;
{\[Alpha], A}=MEObsMinimize[\[Alpha], A] //Chop ;
Return[{\[Alpha], A}];
]


End[(* Private *)];


EndPackage[];
