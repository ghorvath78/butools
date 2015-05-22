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
TestMCPackage::usage = "TestMCPackage[] : Executes various tests to check the functions of the MC package";


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
If[ Not[Q.h-1<=size2 prec],
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


TestMCPackage[]:=Module[{ret,flag,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16,Q17,Q18,Q19,Q20,Q21,Q22,Q23},
Print["---BuTools: MC package test file---"];
Print["Enable the verbose messages with the BuToolsVerbose flag"];
BuTools`Verbose = True;
Print["Enable input parameter checking with the BuToolsCheckInput flag"];
BuTools`CheckInput = True;
On[Assert];
Print["----------------------------"];
?CRPSolve

Print["Input:"];
Print["------"];

Q1= {{-4.3, 3.5, 0.8},{-8.4, 6.5, 1.9},{17.3, -12.7, -4.6}};
Print["Q1=",Q1];

Print["Test:"];
Print["-----"];

Print["CRPSolve[Q1]:"];
ret=CRPSolve[Q1];
Print[ret];
Assert[Norm[ret.Q1]<10^-12, "The solution does not satisfy ret*Q1=0!"];

Print["----------------------------"];
?DRPSolve;

Print["Input:"];
Print["------"];

Q2={{-0.9, 0.5, 1.4},{0.9, -0.9, 1},{0.3, 1.3, -0.6}};
Print["Q2=",Q2];

Print["Test:"];
Print["-----"];

Print["DRPSolve[Q2]:"];
ret=DRPSolve[Q2];
Print[ret];
Assert[Norm[ret.Q2-ret]<10^-12, "The solution does not satisfy ret*Q2=ret!"];

Print["----------------------------"];
?CTMCSolve;

Print["Input:"];
Print["------"];

Q3 = {{0.1, 0.5, 0.4},{0.9, 0.1, 0}};
Print["Q3=",Q3];

Print["Test:"];
Print["-----"];

Print["CTMCSolve[Q3]:"];
Print[Catch[CTMCSolve[Q3]]];

Print["Input:"];
Print["------"];

Q4 = {{0.1, 0.5, 0.4},{0.9, 0.1, 0},{0.9, 0.1, 0}};
Print["Q4=",Q4];

Print["Test:"];
Print["-----"];

Print["CTMCSolve[Q4]:"];
Print[Catch[CTMCSolve[Q4]]];
    
Print["Input:"];
Print["------"];

Q5 = {{-0.9, 0.5, 0.4},{0.9, -0.9, 0},{0.3, 0.3, -0.6}};
Print["Q5=",Q5];

Print["Test:"];
Print["-----"];

Print["CTMCSolve[Q5]:"];
ret=CTMCSolve[Q5];
Print[ret];
Assert[Norm[ret.Q5]<10^-12, "The solution does not satisfy ret*Q5=0!"];

Print["----------------------------"];
?DTMCSolve;

Print["Input:"];
Print["------"];

Q6 = {{0.1, 0.5, 0.4},{0.9, 0.1, 0},{0.3, -0.3, 0.4}};
Print["Q6=",Q6];

Print["Test:"];
Print["-----"];

Print["DTMCSolve[Q6]:"];
Print[Catch[DTMCSolve[Q6]]];

Print["Input:"];
Print["------"];

Q7 = {{0.1, 0.5, 0.4},{0.9, 0.1, 0},{0.3, 0.3, 0.4}};
Print["Q7=",Q7];

Print["Test:"];
Print["-----"];

Print["DTMCSolve[Q7]:"];
ret=DTMCSolve[Q7];
Print[ret];
Assert[Norm[ret.Q7-ret]<10^-12, "The solution does not satisfy ret*Q7=ret!"];

Print["----------------------------"];
?CheckGenerator;

Print["Input:"];
Print["------"];

Q8 = {{-0.9, 0.2, 0.4},{0, 0.9, 0.9},{0, 0.6, -0.6}};

Print["Test:"];
Print["-----"];

Print["CheckGenerator[Q8,true]:"];
flag=CheckGenerator[Q8,True];
Print[flag];
Assert[flag==False,"CheckGenerator did not detect bad row sum!"];

Print["Input:"];
Print["------"];

Q9 = {{-0.9, 0.5, 0.4},{0.9, -0.9, 0},{0.3, 0.3, -0.6}};

Print["Test:"];
Print["-----"];

Print["CheckGenerator[Q9,true]:"];
flag=CheckGenerator[Q9,True];
Print[flag];
Assert[flag==True,"CheckGenerator did not recognize a valid input!"];


Print["Input:"];
Print["------"];

Q10 = {{-0.9, 0.2, 0.4},{0.9, -0.9, 0},{0.3, 0.3, -0.6}};

Print["Test:"];
Print["-----"];

Print["CheckGenerator[Q10,true]:"];
flag=CheckGenerator[Q10,True];
Print[flag];
Assert[flag==True,"CheckGenerator did not recognize a valid input!"];

Print["Input:"];
Print["------"];

Q11 = {{-0.9, 0.5, 0.4},{0.9, -1.1, 0},{0.3, 0.3, -0.6}};

Print["Test:"];
Print["-----"];

Print["CheckGenerator[Q11]:"];
flag=CheckGenerator[Q11];
Print[flag];
Assert[flag==False,"CheckGenerator did not recognize the non-zero row sum!"];

Print["Input:"];
Print["------"];

Q12 = {{-0.9, 0.5, 0.4},{0.9, -0.9, 0},{0.3, 0.3, -0.6}};

Print["Test:"];
Print["-----"];

Print["CheckGenerator[Q12]:"];
flag=CheckGenerator[Q12];
Print[flag];
Assert[flag==True,"CheckGenerator did not recognize a valid input!"];

Print["----------------------------"];
?CheckProbMatrix;

Print["Input:"];
Print["------"];

Q13 = {{0.1, 0.5, 0.4},{0.9, 0.1, 0},{0.3, -0.1, 0.4}};

Print["Test:"];
Print["-----"];

Print["CheckProbMatrix[Q13]:"];
flag=CheckProbMatrix[Q13];
Print[flag];
Assert[flag==False,"CheckProbMatrix did not recognize the negative entry!"];

Print["Input:"];
Print["------"];

Q14 = {{0.1, 0.5, 0.4},{0.9, 0.1, 0},{0.3, 0.1, 0.4}};

Print["Test:"];
Print["-----"];

Print["CheckProbMatrix[Q14]:"];
flag=CheckProbMatrix[Q14];
Print[flag];
Assert[flag==False,"CheckProbMatrix did not recognize the invalid row sum!"];

Print["Input:"];
Print["------"];

Q15 = {{0.1, 0.5, 0.4},{0.9, 0.1, 0},{0.3, 0.3, 0.4}};

Print["Test:"];
Print["-----"];

Print["CheckProbMatrix[Q15]:"];
flag=CheckProbMatrix[Q15];
Print[flag];
Assert[flag==True,"CheckProbMatrix did not recognize that the input is valid!"];

Print["Input:"];
Print["------"];

Q16 = {{0.1, 0.5, 0.4},{0.9, 0.1, 0},{0.3, 0.3, 0.4}};

Print["Test:"];
Print["-----"];

Print["CheckProbMatrix[Q16,true]:"];
flag=CheckProbMatrix[Q16,True];
Print[flag];
Assert[flag==False,"CheckProbMatrix did not recognize wrong transient matrix!"];

Print["Input:"];
Print["------"];

Q17 = {{0.1, 0.5, 0.4},{0.9, 0.1, 0},{0.3, 0.1, 0.4}};

Print["Test:"];
Print["-----"];

Print["CheckProbMatrix[Q17,true]:"];
flag=CheckProbMatrix[Q17,True];
Print[flag];
Assert[flag==True,"CheckProbMatrix did not recognize that the input is valid!"];

Print["----------------------------"];
?CheckProbVector;

Print["Input:"];
Print["------"];

Q18 = {1.1, -0.1};

Print["Test:"];
Print["-----"];

Print["CheckProbVector[Q18]:"];
flag=CheckProbVector[Q18];
Print[flag];
Assert[flag==False,"CheckProbVector did not recognize the negative entry!"];

Print["Input:"];
Print["------"];

Q19 = {1.1, 0.1};

Print["Test:"];
Print["-----"];

Print["CheckProbVector[Q19]:"];
flag=CheckProbVector[Q19];
Print[flag];
Assert[flag==False,"CheckProbVector did not recognize invalid sum!"];

Print["Input:"];
Print["------"];

Q20 = {1, 0};

Print["Test:"];
Print["-----"];

Print["CheckProbVector[Q20]:"];
flag=CheckProbVector[Q20];
Print[flag];
Assert[flag==True,"CheckProbVector did not recognize that the input is valid!"];

Print["Input:"];
Print["------"];

Q21 = {0.9, -0.1};

Print["Test:"];
Print["-----"];

Print["CheckProbVector[Q21,true]:"];
flag=CheckProbVector[Q21,True];
Print[flag];
Assert[flag==False,"CheckProbVector did not recognize the negative entry!"];

Print["Input:"];
Print["------"];

Q22 = {0.9, 0.1};

Print["Test:"];
Print["-----"];

Print["CheckProbVector[Q22,true]:"];
flag=CheckProbVector[Q22,True];
Print[flag];
Assert[flag==True,"CheckProbVector did not recognize that the prob. vector is not transient!"];

Print["Input:"];
Print["------"];

Q23 = {0.8, 0.1}

Print["Test:"];
Print["-----"];

Print["CheckProbVector[Q23,true]:"];
flag=CheckProbVector[Q23,True];
Print[flag];
Assert[flag==True,"CheckProbVector did not recognize that the input is valid!"];
];


End[(* Private *)];
EndPackage[];
