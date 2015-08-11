(* ::Package:: *)

(*
   BuTools RepTrans Package
*)

BeginPackage["BuTools`RepTrans`"];
TransformToAcyclic::usage = "B = TransformToAcyclic[A, maxSize, precision]: Transforms an arbitrary matrix to a Markovian bi-diagonal matrix";
TransformToMonocyclic::usage = "B = TransformToMonocyclic[A, maxSize, precision]: Transforms an arbitrary matrix to a Markovian monocyclic matrix";
ExtendToMarkovian::usage = "{beta, B} = ExtendToMarkovian[alpha, A, maxSize, precision]: Appends an appropriate Erlang tail to the representation that makes the initial vector Markovian";
FindMarkovianRepresentation::usage = "mrep = FindMarkovianRepresentation[rep, transfun, evalfunc, precision]: Obtains a Markovian representation from a non-Markovian one, with keeping the size the same";
SimilarityMatrix::usage = "B = SimilarityMatrix[A1, A2]: Returns the matrix that transforms A1 to A2";
SimilarityMatrixForVectors::usage = "B = SimilarityMatrixForVectors[vecA, vecB]: Returns the similarity transformation matrix that transforms a given column a to a column vector b";
MStaircase::usage = "{B, n} = MStaircase[X, Z, precision]: Computes a smaller representation of RAP using staircase algorithm";


Begin["`Private`"];


If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckInput"]],BuTools`CheckInput=True];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckPrecision"]],BuTools`CheckPrecision=N[10^-12]];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`Verbose"]],BuTools`Verbose=False];


TransformToAcyclic[A_,maxSize_:100,precision_:N[10^-14]] :=
If[Max[Abs[Im[Eigenvalues[A]]]]>=precision,
	Throw["TransformToAcyclic: Complex eigenvalue found, no acyclic representation exists."],
	Return[TransformToMonocyclic[A,maxSize,precision]]
];


TransformToMonocyclic[A_,maxSize_:100,precision_:N[10^-14]] :=
Module[{evalues,febs,size,n,sigma,z,ev,B,pos,Ni},
	(* Generate FEBs *)
	evalues = Eigenvalues[A];
	febs = {};
	size=0;
	Do[
		If[Abs[Im[evi]]<=precision,
			AppendTo[febs,{-Re[evi],0,1,Re[evi]}];
			size++;
		,
			If[Im[evi]<0,Continue[]]; (* complex conjugate pairs are processed at once *)
			(* determine size of the FEB *)
			n=3; size+=3;
			While[-Abs[Im[evi]]/Re[evi]>=Cot[\[Pi]/n],n++;size++; If[size>maxSize,Throw["The represetation is too large (>maxSize). No result returned."]]];
			(* determine rate and feedback prob of the FEB *)
			sigma=-(2 Re[evi]-Abs[Im[evi]](Cot[\[Pi]/n]-Tan[\[Pi]/n]))/2;
			z=(Abs[Im[evi]](Cot[\[Pi]/n]+Tan[\[Pi]/n]) / (2 sigma))^n;
			(* eigenvalues of the FEB. we need the dominant one for sorting *)
			ev=Table[(z^(1/n) Cos[2(k-1)\[Pi]/n]-1)sigma + z^(1/n) sigma Sin[2(k-1)\[Pi]/n]I,{k,n}];
			AppendTo[febs,{sigma,z,n,Max[Re[ev]]}];
		];
	,{evi,evalues}];
	(* Sort FEBs *)
	febs = Sort[febs, #1[[4]]>#2[[4]] &];
	(* Contstruct monocyclic generator *)
	B=ConstantArray[0,{size,size}];
	pos=1;
	Do[
		Ni=feb[[3]];
		B[[pos+Ni-1,pos]]=feb[[2]] feb[[1]];
		Do[B[[pos+i,pos+i]]=-feb[[1]],{i,0,Ni-1}];
		Do[B[[pos+i,pos+i+1]]=feb[[1]],{i,0,Ni-2}];
		If[pos+Ni-1<size,B[[pos+Ni-1,pos+Ni]]=(1-feb[[2]])feb[[1]]];
		pos+=Ni;
	,{feb,febs}];
	Return[B];
];


ExtendToMarkovian[alpha_, A_, maxSize_:100, precision_:N[10^-14]]:=
Module[{t0,t0upper,t0lower,beta,IniVecWithTail,increment,bestT0,bestLupper,success,
	Llower,Lupper,L,B,AN},
	IniVecWithTail[tailLengthx_,mux_]:=
	Module[{len,betax,WG,opv,clv},
		len = Dimensions[A][[1]];
		betax=ConstantArray[0,{len+tailLengthx}];
		WG=IdentityMatrix[len]+A/mux;
		opv=alpha;
		clv=-Total[A/mux,{2}];
		Do[betax[[k]]=opv.clv; opv=opv.WG,{k,len+tailLengthx,len+1,-1}];
		betax[[1;;len]]=opv;
		Return[betax];
	];
	
	(* initial value for t0upper *)
	t0lower=0; t0upper=1;
	beta = alpha.MatrixExp[A t0upper];
	(*instead of Min[beta]<-precision we write AnyTrue[Map[Reduce[Re[#]<-precision]&,beta],TrueQ] since it works with infinite precision*)
	While[AnyTrue[Map[Reduce[Re[#]<-precision]&,beta],TrueQ], t0upper=2*t0upper; beta=alpha.MatrixExp[A t0upper];];
	(* interval bisectioning *)
	While[(t0upper-t0lower)/(t0upper+t0lower)>precision,
		t0=(t0upper+t0lower)/2;
		beta=alpha.MatrixExp[A t0];
		If[AnyTrue[Map[Reduce[Re[#]<-precision]&,beta],TrueQ],t0lower=t0,t0upper=t0];
	];
	t0=t0upper;
    (* find optimal length and rate parameters of the Erlang tail *)
    (* larger t0 sometimes gives fewer states thus we try to increase 
	   t0 gradually till the number of states decreases *)
	increment=11/10;
	bestT0=-1; bestLupper=-1;
	Do[
		Llower=1; Lupper=1;
		beta = IniVecWithTail[Lupper,Lupper/t0];
		While[AnyTrue[Map[Reduce[Re[#]<-precision]&,beta],TrueQ] && Lupper<maxSize, Lupper=2*Lupper; beta = IniVecWithTail[Lupper,Lupper/t0];];
		success=AllTrue[Map[Reduce[Re[#]>=-precision]&,beta],TrueQ];
		If[success,
			While[Lupper-Llower>1,
				L=Round[(Lupper+Llower)/2];
				beta = IniVecWithTail[L,L/t0];
				If[AnyTrue[Map[Reduce[Re[#]<-precision]&,beta],TrueQ],Llower=L,Lupper=L];
			];
		];
		If[success,
			If[bestLupper>=0 && Lupper>bestLupper, Break[], bestLupper=Lupper;bestT0=t0;t0=t0*increment];
		,
			If[bestLupper>=0, Break[], t0=t0*increment];
		];
	,{100}];
	t0=bestT0;
	Lupper=bestLupper;
	If[Lupper<0,Throw["No positive representation found up to the given size!"]];
	(* final result *)
	beta = IniVecWithTail[Lupper,Lupper/t0];
	AN=Dimensions[A][[1]];
	B=ConstantArray[0,{AN+Lupper,AN+Lupper}];
	B[[1;;AN,1;;AN]]=A;
	B[[1;;AN,AN+1]]=-Total[A,{2}];
	Do[
		B[[AN+i,AN+i]]=-Lupper/t0;
		If[i<Lupper,B[[AN+i,AN+i+1]]=Lupper/t0];
	,{i,Lupper}];
	Return[{N[beta,Precision[alpha]],N[B,Precision[A]]}];
];


SimilarityMatrix[A1_,A2_]:=
Module[{N1,N2,Q1,R1,Q2,R2,Q,R,c1,c2,X,M,x,h,Prec,m,sol},
	If[Not[SquareMatrixQ[A1]&&SquareMatrixQ[A2]],Throw["SimilarityMatrix: The input matrices must be square!"]];
	N1=Dimensions[A1][[1]];
	N2=Dimensions[A2][[1]];
	If[N1>N2,Throw["SimilarityMatrix: The first input matrix must be smaller than the second one!"]];
	If[(Precision[A1]==Infinity && Precision[A1]==Infinity) || Not[MatrixQ[A1,NumericQ]] || Not[MatrixQ[A2,NumericQ]],
		X=Table[Subscript[x, i,j],{i,N1},{j,N2}];
		h=ConstantArray[1,{N2,1}];
		sol=Solve[{A1.X==X.A2,X.h==1},Flatten[X]];
		If[Length[sol]==0,Throw["SimilarityMatrix: No solutions found."]];
		Return[X/.sol[[1]]];
	];
	Prec=Min[Precision[A1],Precision[A2]];
	{Q1,R1}=SchurDecomposition[N[A1,Prec],RealBlockDiagonalForm->False];
	{Q2,R2}=SchurDecomposition[N[A2,Prec],RealBlockDiagonalForm->False];
	c1=ConjugateTranspose[{Total[Q2,{1}]}];
	c2=Conjugate[Total[Q1,{1}]];
	X=ConstantArray[0,{N1,N2}];
	Do[
		M=-R2+R1[[k,k]]IdentityMatrix[N2];
		If[k==N1,m=ConstantArray[0,N2],m=-R1[[k,k+1;;]].X[[k+1;;,;;]]];
		{Q,R}=QRDecomposition[Transpose[Join[M,c1,2]]];
		X[[k,;;]]=Inverse[R].Q.Append[m,c2[[k]]];
		(*X[[k,;;]]=LinearSolve[Transpose[Join[M,c1,2]],Append[m,c2[[k]]]];*)
	,{k,N1,1,-1}];
	Return[Re[Q1.X.ConjugateTranspose[Q2]]];
];


SimilarityMatrixForVectors[vecA_,vecB_]:=
Module[{m,ix,P,cpA,B,fVecA,fVecB},
	fVecA=Flatten[vecA];
	fVecB=Flatten[vecB];
	m = Length[fVecA];
	ix=Ordering[-fVecA];
	P=ConstantArray[0,{m,m}];
	Do[P[[i,ix[[i]]]]=1,{i,m}];
	cpA=P.fVecA;
	B=ConstantArray[0,{m,m}];
	Do[B[[i,1;;i]]=fVecB[[i]]/Total[cpA[[1;;i]]],{i,m}];
	Return[B.P];
];


FindMarkovianRepresentation[Rep_, TransFunc_, EvalFunc_, Prec_:N[10^-7]] :=
Module[{NRep, M, k, b, ddist, odist, ErrorMinimize, Elementary},

	Elementary[MRep_, bv_,kv_]:=
	Module[{bestRep,bestDist,repSize,B,newRep,newDist},
		bestDist = EvalFunc[MRep,k];
		bestRep = MRep;
		repSize = Dimensions[Rep[[2]]][[1]];
		Do[
			Do[
				If[i!=j,
					B = IdentityMatrix[repSize];
					B[[i,j]]=b;
					B[[i,i]] = 1-b;
					newRep = TransFunc[MRep,B];
					newDist = EvalFunc[newRep,k];
					If[newDist<bestDist,bestRep=newRep;bestDist=newDist];
					B = IdentityMatrix[repSize];
					B[[i,j]]=-b;
					B[[i,i]] = 1+b;
					newRep = TransFunc[MRep,B];
					newDist = EvalFunc[newRep,k];
					If[newDist<bestDist,bestRep=newRep;bestDist=newDist];
				];
			,{j,repSize}];
		,{i,repSize}];
		Return[{bestRep,bestDist}];
	];

	ErrorMinimize[MRep_,bv_,kv_]:=
	Module[{bestRep,currRep,dist,lastDist},
		lastDist = EvalFunc[MRep,k];
		bestRep = MRep;
		currRep = MRep;
		Do[
			{currRep, dist}=Elementary[currRep, bv, kv];
			If[dist>=lastDist,Break[],lastDist=dist;bestRep=currRep];
		,{i,M^2}];
		Return[{bestRep,lastDist}];
	];

	NRep=Rep;
	M=Dimensions[Rep[[2]]][[1]];
	b=N[1/2,Precision[Rep]];
	odist=Infinity;
	While[b>Prec/2,
		Do[
			Do[
				{NRep,ddist}=ErrorMinimize[NRep,b,k];
				If[ddist<Prec, Return[NRep];];
			,{k,4}];
			If[odist<=ddist,Break[],odist=ddist];
		,{m,M^2}];
		b=b/2;
	];
	Return[NRep];
];


MStaircase[Xin_,Zin_,prec_: N[10^-8]]:=
Module[{X,Z,msize,m,U,ranksum,crit,r,Ui,S,T,TEMP,Z1,X1,X2,X3,X4,Transf,I,n,x,i,nonzero,zeros,R,
y,\[CapitalGamma],TEMP1,sizeUi,mm,ii},

X=Xin;Z=Zin;
msize=Length[X];
m=Dimensions[X[[1]]][[1]];
U=IdentityMatrix[m];
ranksum=0;
crit=True;
While[crit, 
	r=MatrixRank[Z,Tolerance->prec];
	ranksum+=r;
	{Ui,S,T}=SingularValueDecomposition[Z];
	(*Calculation of the new U,X,Z matrices*)
	sizeUi=Dimensions[Ui][[1]];
	Transf=IdentityMatrix[ranksum-r+sizeUi];
	Transf[[-sizeUi ;;-1,-sizeUi ;;-1]]=Transpose[Ui];
	U=Transpose[Transf.Transpose[U]];

	Z=ConstantArray[0,{Dimensions[Ui][[1]]-r,0}];
	Do[
		TEMP=Transpose[Ui].X[[ii]].Ui;
		X[[ii]]=TEMP[[r+1;;-1,r+1;;-1]];
		Z=Join[Z,TEMP[[r+1;;-1,1;;r]],2];
	,{ii,msize}];
	If[Norm[Z]<prec || MatrixRank[Z,Tolerance->prec]==m-ranksum, crit=False];
];

n=ranksum;
I=IdentityMatrix[m];
If[Norm[Z]<prec,
	n=ranksum;
	x=Total[Transpose[U],{2}];
	x=x[[1;;n]];

	R=IdentityMatrix[n];
	If[Min[Abs[x]]<prec,
		nonzero=0;
		zeros=Table[False,{n}];
		Do[		
			If[Abs[x[[i]]]<prec,zeros[[i]]=True,If[nonzero==0,nonzero=i]];
		,{i,n}];
		If[nonzero!=0,
			Do[If[zeros[[i]],R[[i,nonzero]]=1],{i,n}];
		]; 
	];
	y=R.x;
	\[CapitalGamma]=DiagonalMatrix[y];
	TEMP=I;
	TEMP[[1;;n,1;;n]]=Inverse[\[CapitalGamma]];
	TEMP1=I;
	TEMP1[[1;;n,1;;n]]=R;
	Return[{Inverse[TEMP.TEMP1.Transpose[U]],n}];
];
If[MatrixRank[Z,Tolerance->prec]==m-ranksum,Return[{I,m}]];
];


End[(* Private *)];
EndPackage[];
