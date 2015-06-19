(* ::Package:: *)

(*
   BuTools MC Package
*)

BeginPackage["BuTools`MC`"];
CTMCSolve::usage = "pi = CTMCSolve[Q]: Returns the steady state solution of a continuous time Markov chain";
DTMCSolve::usage = "pi = DTMCSolve[Q]: Returns the steady state solution of a discrete time Markov chain";
CRPSolve::usage = "pi = CRPSolve[Q]: Returns the steady state solution of a continuous time rational process";
DRPSolve::usage = "pi = DRPSolve[Q]: Returns the steady state solution of a discrete time rational process.";
CheckGenerator::usage = "r = CheckGenerator[Q, transient, prec]: Checks if a matrix is a valid generator of a CTMC";
CheckProbMatrix::usage = "r = CheckProbMatrix[P, transient, prec]: Checks if a matrix is a valid transition probability matrix of a DTMC";
CheckProbVector::usage = "r = CheckProbVector[pi, sub, prec]: Checks if a vector is a valid probability vector";


Begin["`Private`"];


If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckInput"]],BuTools`CheckInput=True];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckPrecision"]],BuTools`CheckPrecision=N[10^-12]];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`Verbose"]],BuTools`Verbose=False];


CRPSolve[Q_] :=
Module[{mx},
	If[BuTools`CheckInput && Max[Abs[Total[Q,{2}]]]>BuTools`CheckPrecision,Throw["CRPSolve: The matrix has a rowsum which isn't zero!"]];
	mx=Q;
	mx[[;;,1]]=1;
	Return[LinearSolve[Transpose[mx], Prepend[ConstantArray[0,Dimensions[Q][[1]]-1],1]]];
]


DRPSolve[P_] :=(
If[BuTools`CheckInput && Max[Abs[Total[P,{2}]]-1]>BuTools`CheckPrecision,Throw["DRPSolve: The matrix has a rowsum which isn't 1!"]];
Return[CRPSolve[P-IdentityMatrix[Dimensions[P][[1]]]]];)


CTMCSolve[Q_] :=
(If[BuTools`CheckInput && Not[CheckGenerator[Q,False]],Throw["CTMCSolve/DTMCSolve: The given matrix is not a valid generator. If you are sure you want this, use CRPSolve/DRPSolve instead CTMCSolve/DTMCSolve."]];
Return[CRPSolve[Q]]; )


DTMCSolve[P_]:= CTMCSolve[P-IdentityMatrix[Dimensions[P][[1]]]];


CheckGenerator[Q_,transient_ : False,precision_ : Null]:=
Module[ {k,myh0,h,size1,size2,prec},
If[Not[NumericQ[precision]], prec=BuTools`CheckPrecision, prec=precision];
{size1,size2}=Dimensions[Q];
If[size1!=size2,
If[BuTools`Verbose==True,Print["CheckGenerator: the generator is not a square matrix!"]];
Return[False]
];

myh0=Q;
For[k=1,k<=size1,k++,
If[myh0[[k,k]]>=prec,
If[BuTools`Verbose==True,Print["CheckGenerator: The diagonal of the generator is not negative (at precision ",prec,")!"]];
Return[False],
myh0[[k,k]]=0
];
];

If[ Min[myh0]<-prec,
If[BuTools`Verbose==True,Print["CheckGenerator: the generator has negative off-diagonal element (at precision ",prec,")!"]];
Return[False]
];

h=Table[1,{size1}];
If[transient,
If[Max[Q.h]>prec,
If[BuTools`Verbose==True,Print["CheckGenerator: A rowsum of the transient generator is greather than 0 (at precision ",prec,")!"]];
Return[False];
];
If[Max[Re[Eigenvalues[Q//N]]]>=prec,
If[BuTools`Verbose==True,Print["CheckGenerator: the transient generator has non-negative eigenvalue (at precision ",prec,")!"]];
Return[False]
],
If[  Max[Abs[Q.h]]>prec, 
If[BuTools`Verbose==True,Print["CheckGenerator: A rowsum of the generator is not 0 (precision:",prec,")!!"]];
Return[False]];
];

Return[True]
];


CheckProbMatrix[Q_,transient_ : False,precision_ :Null]:=
Module[{size1, size2, h, prec},
If[Not[NumericQ[precision]], prec=BuTools`CheckPrecision, prec=precision];
{size1,size2}=Dimensions[Q];
If[size1!=size2,
If[BuTools`Verbose==True,Print["CheckProbMatrix: the matrix is not a square matrix!"]];
Return[False]
];

If[ Min[Q]<-prec,
If[BuTools`Verbose==True,Print["CheckProbMatrix: the matrix has negative element (at precision ",prec,")!"]];
Return[False]
];

h=Table[1,{size1}];

If[transient,
If[ Max[Q.h-h]>size2 prec,
If[BuTools`Verbose==True,
Print["CheckProbMatrix: A rowsum of the transient matrix is not less or equal than 1!"]];
Return[False]
];
If[Not[Max[Re[Eigenvalues[Q//N]]]<1-prec],
If[BuTools`Verbose==True,
Print["CheckProbMatrix: The real part of the largest eigenvalue of the transient matrix is not less than 1!"]];
Return[False]
],
If[ Max[Abs[Q.h-h]]>size2 prec,
If[BuTools`Verbose==True,Print["CheckProbMatrix: A rowsum of the matrix is not 1 (precision:",prec,")!!"]]; 
Return[False]
]
];

Return[True]
];


CheckProbVector[pi_,sub_ : False, precision_ : Null]:=
Module[{prec},
If[Not[NumericQ[precision]], prec=BuTools`CheckPrecision, prec=precision];
If[ Min[pi]<-prec,
If[BuTools`Verbose==True,Print["CheckProbVector: The vector has negative element!"]];
Return[False]
];

If[sub,
If[ Not[Total[pi]<1+prec Length[pi]],
If[BuTools`Verbose==True,
Print["CheckProbVector: The sum of the substochastic vector is not less than 1!"];
];
Return[False]
],
If[ Abs[Total[pi]-1]>prec Length[pi],
If[BuTools`Verbose==True,Print["CheckProbVector: The sum of the vector is not 1 (precision:",prec,")!"]];
Return[False]
]
];
Return[True]
];


End[(* Private *)];
EndPackage[];
