(* ::Package:: *)

(*
   BuTools Queues Package
*)

BeginPackage["BuTools`Queues`"];
QBDQueue::usage="";
MAPMAP1::usage="";
FluidQueue::usage="";
FluFluQueue::usage="";
MMAPPH1PRPR::usage="";
MMAPPH1NPPR::usage="";


Begin["`Private`"];


Needs["BuTools`MC`"];
Needs["BuTools`MAM`"];
Needs["BuTools`RepTrans`"];
Needs["BuTools`Moments`"];
Needs["BuTools`MAP`"];
Needs["BuTools`PH`"];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckInput"]],BuTools`CheckInput=True];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckPrecision"]],BuTools`CheckPrecision=N[10^-12]];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`Verbose"]],BuTools`Verbose=False];


QBDQueue[B_,L_,F_,L0_,argv__]:=
Module[{varargin,prec,needST,eaten,
pi0,R,U,Rh,NN,II,eta,z,Ret,retIx,argIx,alpha,A,Bi,numOfMoms,moms,
Z,iZ,points,Bm,Bmi,iR,Delta,ix,nz},
	varargin = List[argv];
    (* parse options *)
    prec = N[10^-14];
    needST = 0;
    eaten = {};
	Do[
		If[varargin[[i]]=="prec",
			prec = varargin[[i+1]];
			AppendTo[eaten,i]; 
			AppendTo[eaten,i+1];
		, If[StringQ[varargin[[i]]] && StringLength[varargin[[i]]]>2 && StringTake[varargin[[i]],2]=="st",
			needST = True;
		]]
	,{i,Length[varargin]}];

    If[BuTools`CheckInput && Not[CheckGenerator[B+L+F]], Throw["QBDQueue: The matrix sum (B+L+F) is not a valid generator of a Markov chain!"]]; 
    If[BuTools`CheckInput && Not[CheckGenerator[L0+F]], Throw["QBDQueue: The matrix sum (L0+F) is not a valid generator of a Markov chain!"]];

    {pi0, R} = QBDSolve[B, L, F, L0, prec];
    NN = Length[pi0];
    II = IdentityMatrix[NN];
    
    If[needST,
        U = L + R.B;
        Rh = Inverse[-U].F;
        eta = pi0.F.Inverse[II-Rh];
        eta = eta/Total[eta];
        z = ArrayReshape[II,{NN*NN,1}];
    ];
    Ret = {};
    argIx = 1;
    While[argIx<=Length[varargin],
        If[MemberQ[eaten, argIx],
            argIx++;
            Continue[];
		,If[varargin[[argIx]]=="qlDistrDPH",
            (* transform it to DPH *)
            alpha = pi0.R.Inverse[II-R];
            A = Inverse[DiagonalMatrix[alpha]].Transpose[R].DiagonalMatrix[alpha];
            AppendTo[Ret,{alpha, A}];
        ,If[varargin[[argIx]]=="qlDistrMG",
            (* transform it to MG *)
            Bm = SimilarityMatrixForVectors[Total[Inverse[II-R].R,{2}], Table[1,{NN},{1}]];
            Bi = Inverse[Bm];
            A = Bm.R.Bi;
            alpha = pi0.Bi;       
            AppendTo[Ret,{alpha, A}];
        ,If[varargin[[argIx]]=="qlMoms",
            numOfMoms = varargin[[argIx+1]];
            argIx++;
            iR = Inverse[II-R];
			moms = Table[m! Total[pi0.MatrixPower[iR,m+1].MatrixPower[R,m]],{m,numOfMoms}];
            AppendTo[Ret,MomsFromFactorialMoms[moms]];
        ,If[varargin[[argIx]]=="qlDistr",
            points = varargin[[argIx+1]];
            argIx++;
			AppendTo[Ret,Table[Total[pi0.If[po==0,IdentityMatrix[NN],MatrixPower[R,po]]],{po,points}]];
        ,If[varargin[[argIx]]=="stDistrPH",
            (* transform to ph distribution *)
            nz = Flatten[Position[eta,_?(#>prec&)]];
            Delta = DiagonalMatrix[eta];
            A = KroneckerProduct[L+F,II[[nz,nz]]] + KroneckerProduct[B,Inverse[Delta[[nz,nz]]].Transpose[Rh[[nz,nz]]].Delta[[nz,nz]]];
            alpha = Flatten[z].KroneckerProduct[II,Delta[[All,nz]]];
            AppendTo[Ret,{alpha, A}];
        ,If[varargin[[argIx]]=="stDistrME",
            (* transform it such that the closing vector is a vector of ones *)
            (* this is the way butools accepts ME distributions *)
            Bm = SimilarityMatrixForVectors[z,Table[1,{Length[z]},{1}]];
            Bmi = Inverse[Bm];
            A = Bm.(KroneckerProduct[Transpose[L+F],II] + KroneckerProduct[Transpose[B],Rh]).Bmi;
            alpha = (KroneckerProduct[Table[1,{1},{NN}], eta].Bmi)[[1]];
            AppendTo[Ret,{alpha, A}];
        ,If[varargin[[argIx]]=="stMoms",
            numOfMoms = varargin[[argIx+1]];
            argIx++;
            Z = KroneckerProduct[Transpose[L+F],II] + KroneckerProduct[Transpose[B],Rh];
            iZ = Inverse[-Z];
			AppendTo[Ret,Table[m!*Total[KroneckerProduct[Table[1,{1},{NN}], eta].MatrixPower[iZ,m+1].(-Z).z,2],{m,numOfMoms}]];
        ,If[varargin[[argIx]]=="stDistr",
            points = varargin[[argIx+1]];
            argIx++;
            Z = KroneckerProduct[Transpose[L+F],II] + KroneckerProduct[Transpose[B],Rh];
            AppendTo[Ret,Table[1-Total[KroneckerProduct[Table[1,{1},{NN}], eta].MatrixExp[Z*po].z,2],{po,points}]];
        ,
            Throw[StringJoin["QBDQueue: Unknown parameter ", ToString[varargin[[argIx]]]]];
        ]]]]]]]]];
        argIx++;
    ];
    If[Length[Ret]==1, Return[Ret[[1]]],Return[Ret]];
];


MAPMAP1[D0_,D1_,S0_,S1_,argv__]:=
Module[{varargin,prec,eaten,needST,IA,IS,B,L,F,L0,pi0,R,NN,II,
U,T,eta,Rh,argIx,Ret,Bm,Bi,numOfMoms,moms,iR,points,theta,nz,vv,
Delta,iT,alpha,A,beta},
	varargin = List[argv];
    (* parse options *)
    prec = N[10^-14];
    needST = 0;
    eaten = {};
	Do[
		If[varargin[[i]]=="prec",
			prec = varargin[[i+1]];
			AppendTo[eaten,i]; 
			AppendTo[eaten,i+1];
		, If[StringQ[varargin[[i]]] && StringLength[varargin[[i]]]>2 && StringTake[varargin[[i]],2]=="st",
			needST = True;
		]]
	,{i,Length[varargin]}];

    If[BuTools`CheckInput && Not[CheckMAPRepresentation[D0,D1]], Throw["MAPMAP1: The arrival process (D0,D1) is not a valid MAP representation!"]]; 
    If[BuTools`CheckInput && Not[CheckMAPRepresentation[S0,S1]], Throw["MAPMAP1: The service process (S0,S1) is not a valid MAP representation!"]];

    IA = IdentityMatrix[Length[D0]];
    IS = IdentityMatrix[Length[S0]];
    
    B = KroneckerProduct[IA,S1];
    L = KroneckerProduct[D0,IS]+KroneckerProduct[IA,S0];
    F = KroneckerProduct[D1,IS];
    L0 = KroneckerProduct[D0,IS];

    {pi0, R} = QBDSolve[B, L, F, L0, prec];
    NN = Length[pi0];
    II = IdentityMatrix[NN];
    
    If[needST,
        (* calculate the distribution of the age at departures *)
        U = L + R.B;
        Rh = Inverse[-U].F;
        T = KroneckerProduct[IA,S0] + Rh.B;
        eta = pi0.F.Inverse[II-Rh];
        eta = eta/Total[eta];
	];

    Ret = {};
    argIx = 1;
    While[argIx<=Length[varargin],
		If[MemberQ[eaten, argIx],
            argIx++;
            Continue[];
		,If[varargin[[argIx]]=="qlDistrDPH",
            (* transform it to DPH *)
            alpha = pi0.R.Inverse[II-R];
            A = Inverse[DiagonalMatrix[alpha]].Transpose[R].DiagonalMatrix[alpha];
            AppendTo[Ret,{alpha, A}];
        ,If[varargin[[argIx]]=="qlDistrMG",
            (* transform it to MG *)
            Bm = SimilarityMatrixForVectors[Total[Inverse[II-R].R,{2}], Table[1,{NN},{1}]];
            Bi = Inverse[Bm];
            A = Bm.R.Bi;
            alpha = pi0.Bi;       
            AppendTo[Ret,{alpha, A}];
        ,If[varargin[[argIx]]=="qlMoms",
            numOfMoms = varargin[[argIx+1]];
            argIx++;
            iR = Inverse[II-R];
			moms = Table[m! Total[pi0.MatrixPower[iR,m+1].MatrixPower[R,m]],{m,numOfMoms}];
            AppendTo[Ret,MomsFromFactorialMoms[moms]];
        ,If[varargin[[argIx]]=="qlDistr",
            points = varargin[[argIx+1]];
            argIx++;
			AppendTo[Ret,Table[Total[pi0.If[po==0,IdentityMatrix[NN],MatrixPower[R,po]]],{po,points}]];
        ,If[varargin[[argIx]]=="stDistrPH",
            (* transform it to PH representation *)
            beta = CTMCSolve[S0+S1];
            theta = DTMCSolve[Inverse[-D0].D1];
            vv = KroneckerProduct[{theta},{beta}][[1]];
            nz = Flatten[Position[vv,_?(#>prec&)]];
            Delta = DiagonalMatrix[vv[[nz]]];
            alpha = Table[1,{NN}].Transpose[B[[nz,All]]].Delta / Total[beta.S1];
            A = Inverse[Delta].Transpose[T[[nz,nz]]].Delta;
            AppendTo[Ret,{alpha, A}];        
        ,If[varargin[[argIx]]=="stDistrME",
            AppendTo[Ret,{eta, T}];
        ,If[varargin[[argIx]]=="stMoms",
            numOfMoms = varargin[[argIx+1]];
            argIx++;
			iT = Inverse[-T];
            AppendTo[Ret,Table[m!*Total[eta.MatrixPower[iT,m]],{m,numOfMoms}]];
        ,If[varargin[[argIx]]=="stDistr",
            points = varargin[[argIx+1]];
            argIx++;
            AppendTo[Ret,Table[1-Total[eta.MatrixExp[T*po]],{po,points}]];
        ,
            Throw[StringJoin["MAPMAP1: Unknown parameter ", ToString[varargin[[argIx]]]]];
        ]]]]]]]]];
        argIx++;
    ];
    If[Length[Ret]==1, Return[Ret[[1]]],Return[Ret]];
];


FluidQueue[Q_,Rin_,Rout_,argv__]:=
Module[{varargin,prec,needST,eaten,Q0,mass0,ini,K,clo,NN,Qm,Rm,
iniKi,lambda,Ret,argIx,Delta,A,alpha,Bm,Bi,numOfMoms,iK,points,
Z,iZ,kini,kclo},
	varargin = List[argv];
    (* parse options *)
    prec = N[10^-14];
    needST = 0;
	Q0 = {};
    eaten = {};
	Do[
		If[varargin[[i]]=="prec",
			prec = varargin[[i+1]];
			AppendTo[eaten,i]; 
			AppendTo[eaten,i+1];
		, If[varargin[[i]]=="Q0",
			Q0 = varargin[[i+1]];
			AppendTo[eaten,i]; 
			AppendTo[eaten,i+1];
		, If[StringQ[varargin[[i]]] && StringLength[varargin[[i]]]>2 && StringTake[varargin[[i]],2]=="st",
			needST = True;
		]]]
	,{i,Length[varargin]}];

    If[BuTools`CheckInput && Not[CheckGenerator[Q,False]], Throw["FluidQueue: Generator matrix Q is not Markovian!"]]; 
    If[BuTools`CheckInput && Length[Q0]>0 && Not[CheckGenerator[Q0]], Throw["FluidQueue: Generator matrix Q0 is not Markovian!"]];
	If[BuTools`CheckInput && (AnyTrue[Diagonal[Rin],#<-BuTools`CheckPrecision &] || AnyTrue[Diagonal[Rout],#<-BuTools`CheckPrecision &]),Throw["FluidQueue: Fluid rates Rin and Rout must be non-negative !"]];

    {mass0, ini, K, clo} = GeneralFluidSolve[Q, Rin-Rout, Q0, prec];
    If[needST,
        NN = Dimensions[Q][[1]];
		{Qm,Rm}=QRDecomposition[Transpose[K]];
		iniKi=Flatten[Inverse[Rm].Qm.(-ini)];
        lambda = Total[mass0.Rin + iniKi.clo.Rin];
    ];
    Ret = {};
    argIx = 1;
    While[argIx<=Length[varargin],
        If[MemberQ[eaten, argIx],
            argIx++;
            Continue[];
		,If[varargin[[argIx]]=="qlDistrPH",
            (* transform it to PH *)
			{Qm,Rm}=QRDecomposition[Transpose[K]];
			iniKi=Flatten[Inverse[Rm].Qm.(-ini)];
            Delta = DiagonalMatrix[iniKi]; (* Delta = diag (ini*inv(-K)); *)
            A = Inverse[Delta].Transpose[K].Delta;
            alpha = Total[clo,{2}].Delta;
            AppendTo[Ret,{alpha, A}];
        ,If[varargin[[argIx]]=="qlDistrME",
            (* transform it to ME *)
            Bm = SimilarityMatrixForVectors[Inverse[-K].Total[clo,{2}], Table[1,{Length[K]},{1}]];
            Bi = Inverse[Bm];
            alpha = ini.Bi;
            A = Bm.K.Bi;
			AppendTo[Ret,{alpha, A}];
        ,If[varargin[[argIx]]=="qlMoms",
            numOfMoms = varargin[[argIx+1]];
            argIx++;
            iK = Inverse[-K];
			AppendTo[Ret,Table[m! Total[ini.MatrixPower[iK,m+1].clo],{m,numOfMoms}]];
        ,If[varargin[[argIx]]=="qlDistr",
            points = varargin[[argIx+1]];
            argIx++;
            iK = Inverse[-K];
			AppendTo[Ret,Table[Total[mass0]+Total[ini.(IdentityMatrix[Length[K]]-MatrixExp[K*po]).iK.clo],{po,points}]];
        ,If[varargin[[argIx]]=="stDistrPH",
            (* transform to ph distribution *)
            Delta = DiagonalMatrix[iniKi/lambda];
            alpha = Flatten[Transpose[clo.Rin]].KroneckerProduct[IdentityMatrix[NN],Delta];
            A = KroneckerProduct[Rout, Inverse[Delta].Transpose[K].Delta] + KroneckerProduct[Q, IdentityMatrix[Length[K]]];
            AppendTo[Ret,{alpha, A}];
        ,If[varargin[[argIx]]=="stDistrME",
            (* transform it such that the closing vector is a vector of ones *)
            (* this is the way butools accepts ME distributions *)
            Bm = SimilarityMatrixForVectors[Transpose[ArrayReshape[Transpose[Inverse[-K].clo.Rin],{1,NN*Length[ini]}]],Table[1,{NN*Length[ini]},{1}]];
            Bi = Inverse[Bm];
            alpha = KroneckerProduct[Table[1,{1},{NN}], {ini/lambda}][[1]].Bi;
            A = Bm.(KroneckerProduct[Transpose[Q],IdentityMatrix[Length[K]]] + KroneckerProduct[Rout,K]).Bi;        
            AppendTo[Ret,{alpha, A}];
        ,If[varargin[[argIx]]=="stMoms",
            numOfMoms = varargin[[argIx+1]];
            argIx++;
			Z = KroneckerProduct[Transpose[Q],IdentityMatrix[Length[K]]] + KroneckerProduct[Rout,K];
            iZ = Inverse[-Z];
            kini = KroneckerProduct[Table[1,{1},{NN}], {ini/lambda}][[1]];
            kclo = Transpose[ArrayReshape[Transpose[Inverse[-K].clo.Rin],{1,NN*Length[ini]}]];
			AppendTo[Ret,Table[m!*Total[kini.MatrixPower[iZ,m+1].(-Z).kclo,2],{m,numOfMoms}]];
        ,If[varargin[[argIx]]=="stDistr",
            points = varargin[[argIx+1]];
            argIx++;
			Z = KroneckerProduct[Transpose[Q],IdentityMatrix[Length[K]]] + KroneckerProduct[Rout,K];
            iZ = Inverse[-Z];
            kini = KroneckerProduct[Table[1,{1},{NN}], {ini/lambda}][[1]];
            kclo = Transpose[ArrayReshape[Transpose[Inverse[-K].clo.Rin],{1,NN*Length[ini]}]];
            AppendTo[Ret,Table[1-Total[kini.MatrixExp[Z*po].kclo,2],{po,points}]];
        ,
            Throw[StringJoin["FluidQueue: Unknown parameter ", ToString[varargin[[argIx]]]]];
        ]]]]]]]]];
        argIx++;
    ];
    If[Length[Ret]==1, Return[Ret[[1]]],Return[Ret]];
];


FluFluQueue[Qin_,Rin_,Qout_,Rout_,srv0stop_,argv__]:=
Module[{varargin,prec,needST,needQL,eaten,Iin,Iout,Q,Q0,mass0,ini,K,clo,NN,Qm,Rm,
iniKi,lambda,mu,Ret,argIx,Delta,A,alpha,Bm,Bi,numOfMoms,iK,points,
kclo,Rh,Qh,massh,inih,Kh,iKh,cloh},
	varargin = List[argv];
    (* parse options *)
    prec = N[10^-14];
    needST = 0;
	needQL = 0;
    eaten = {};
	Do[
		If[varargin[[i]]=="prec",
			prec = varargin[[i+1]];
			AppendTo[eaten,i]; 
			AppendTo[eaten,i+1];
		, If[StringQ[varargin[[i]]] && StringLength[varargin[[i]]]>2 && StringTake[varargin[[i]],2]=="st",
			needST = True;
		, If[StringQ[varargin[[i]]] && StringLength[varargin[[i]]]>2 && StringTake[varargin[[i]],2]=="ql",
			needQL = True;
		]]]
	,{i,Length[varargin]}];

    If[BuTools`CheckInput && Not[CheckGenerator[Qin,False]], Throw["FluFlu: Generator matrix Qin is not Markovian!"]]; 
    If[BuTools`CheckInput && Not[CheckGenerator[Qout,False]], Throw["FluFlu: Generator matrix Qout is not Markovian!"]];
	If[BuTools`CheckInput && (AnyTrue[Diagonal[Rin],#<-BuTools`CheckPrecision &] || AnyTrue[Diagonal[Rout],#<-BuTools`CheckPrecision &]),Throw["FluFluQueue: Fluid rates Rin and Rout must be non-negative !"]];

	Iin = IdentityMatrix[Length[Qin]];
	Iout = IdentityMatrix[Length[Qout]];
	If[needQL,
		Q = KroneckerProduct[Qin,Iout] + KroneckerProduct[Iin,Qout];
		If[srv0stop,
			Q0 = KroneckerProduct[Qin,Iout] + KroneckerProduct[Rin,PseudoInverse[Rout].Qout];
		,
			Q0 = {};
		];
	    {mass0, ini, K, clo} = GeneralFluidSolve[Q, KroneckerProduct[Rin,Iout]-KroneckerProduct[Iin,Rout], Q0, prec];
	];
    If[needST,
        Rh = KroneckerProduct[Rin,Iout] - KroneckerProduct[Iin,Rout];
        Qh = KroneckerProduct[Qin, Rout] + KroneckerProduct[Rin, Qout];       
        {massh, inih, Kh, cloh} = GeneralFluidSolve[Qh, Rh, {}, prec];
        (* sojourn time density in case of *)
        (* srv0stop = false: inih*expm(Kh*x)*cloh*kron(Rin,Iout)/lambda *)
        (* srv0stop = true: inih*expm(Kh*x)*cloh*kron(Rin,Rout)/lambda/mu *)
        lambda = Total[CTMCSolve[Qin].Rin];
        mu = Total[CTMCSolve[Qout].Rout];
		If[Not[srv0stop],
			kclo = cloh.KroneckerProduct[Rin,Iout]/lambda;
		,
			kclo = cloh.KroneckerProduct[Rin,Rout]/lambda/mu;
		];
    ];
    Ret = {};
    argIx = 1;
    While[argIx<=Length[varargin],
        If[MemberQ[eaten, argIx],
            argIx++;
            Continue[];
		,If[varargin[[argIx]]=="qlDistrPH",
            (* transform it to PH *)
			{Qm,Rm}=QRDecomposition[Transpose[K]];
			iniKi=Flatten[Inverse[Rm].Qm.(-ini)];
            Delta = DiagonalMatrix[iniKi]; (* Delta = diag (ini*inv(-K)); *)
            A = Inverse[Delta].Transpose[K].Delta;
            alpha = Total[clo,{2}].Delta;
            AppendTo[Ret,{alpha, A}];
        ,If[varargin[[argIx]]=="qlDistrME",
            (* transform it to ME *)
            Bm = SimilarityMatrixForVectors[Inverse[-K].Total[clo,{2}], Table[1,{Length[K]},{1}]];
            Bi = Inverse[Bm];
            alpha = ini.Bi;
            A = Bm.K.Bi;
			AppendTo[Ret,{alpha, A}];
        ,If[varargin[[argIx]]=="qlMoms",
            numOfMoms = varargin[[argIx+1]];
            argIx++;
            iK = Inverse[-K];
			AppendTo[Ret,Table[m! Total[ini.MatrixPower[iK,m+1].clo],{m,numOfMoms}]];
        ,If[varargin[[argIx]]=="qlDistr",
            points = varargin[[argIx+1]];
            argIx++;
            iK = Inverse[-K];
			AppendTo[Ret,Table[Total[mass0]+Total[ini.(IdentityMatrix[Length[K]]-MatrixExp[K*po]).iK.clo],{po,points}]];
        ,If[varargin[[argIx]]=="stDistrPH",
            (* transform to ph distribution *)
			{Qm,Rm}=QRDecomposition[Transpose[Kh]];
			iniKi=Flatten[Inverse[Rm].Qm.(-inih)];
            Delta = DiagonalMatrix[iniKi];
			alpha = Total[Delta.kclo,{2}];
			A = Inverse[Delta].Transpose[Kh].Delta;
            AppendTo[Ret,{alpha, A}];
        ,If[varargin[[argIx]]=="stDistrME",
            (* transform it such that the closing vector is a vector of ones *)
            (* this is the way butools accepts ME distributions *)
			Bm = SimilarityMatrixForVectors[Transpose[{Total[kclo,{2}]}],Table[1,{Length[Kh]},{1}]];
			Bi = Inverse[Bm];
            alpha = inih.Inverse[-Kh].Bi;
            A = Bm.Kh.Bi;        
            AppendTo[Ret,{alpha, A}];
        ,If[varargin[[argIx]]=="stMoms",
            numOfMoms = varargin[[argIx+1]];
            argIx++;
			iKh = Inverse[-Kh];
            AppendTo[Ret,Table[m!*Total[inih.MatrixPower[iKh,m+1].kclo],{m,numOfMoms}]];
        ,If[varargin[[argIx]]=="stDistr",
            points = varargin[[argIx+1]];
            argIx++;
			iKh = Inverse[-Kh];
            AppendTo[Ret,Table[1-Total[inih.MatrixExp[Kh*po].iKh.kclo],{po,points}]];
        ,
            Throw[StringJoin["FluFluQueue: Unknown parameter ", ToString[varargin[[argIx]]]]];
        ]]]]]]]]];
        argIx++;
    ];
    If[Length[Ret]==1, Return[Ret[[1]]],Return[Ret]];
];


MMAPPH1PRPR[Dm_,sigma_,S_,argv__]:=
Module[{varargin,prec,eaten,classes,K,erlMaxOrder,
D0,ALMA,s,sD,M,NN,II,Ret,sM,Qwmm,Qwmp,Qwpm,Qwpp,kix,
Psiw,Kw,Uw,Ua,b,X,pm,Qm,Rm,Bw,kappa,argIx,
numOfSTMoms,res,A,B,C,bino,Qsmm,Qsmp,Qspm,Qspp,Np,
inis,Psis,P,Pn,rtMoms,stCdfPoints,pr,lambda,L,Psie,
numOfQLMoms,QLDPn,dqlMoms,lambdak,pi,QLPn,iTerm,sumP,
qlMoms,numOfQLProbs,sDk,Psid,dqlProbs,qlProbs,retIx,
F,L0,p0,R,Rpow},
	varargin = List[argv];
	K = Length[Dm]-1;
    (* parse options *)
    prec = N[10^-14];
    erlMaxOrder = 200;
	classes = Range[K];
    eaten = {};
	Do[
		If[varargin[[i]]=="prec",
			prec = varargin[[i+1]];
			AppendTo[eaten,i]; 
			AppendTo[eaten,i+1];
		, If[varargin[[i]]=="erlMaxOrder",
			erlMaxOrder = varargin[[i+1]];
			AppendTo[eaten,i]; 
			AppendTo[eaten,i+1];
		, If[varargin[[i]]=="classes",
			classes = varargin[[i+1]];
			AppendTo[eaten,i]; 
			AppendTo[eaten,i+1];
		]]];
	,{i,Length[varargin]}];

    If[BuTools`CheckInput && Not[CheckMMAPRepresentation[Dm]], Throw["MMAPPH1PRPR: The arrival process is not a valid MMAP representation!"]]; 
    Do[If[BuTools`CheckInput && Not[CheckPHRepresentation[sigma[[k]],S[[k]]]], Throw["MMAPPH1PRPR: the vector and matrix describing the service times is not a valid PH representation!"]],{k,K}];

    (* some preparation *)
    D0 = Dm[[1]];
    NN = Dimensions[D0][[1]];
    II = IdentityMatrix[NN];
    sD = Sum[Dm[[i]],{i,K+1}];
	s = Table[Transpose[{Total[-S[[i]],{2}]}],{i,K}];
	M = Table[Length[sigma[[i]]],{i,K}];
    Ret = {};
	Do[
        (* step 1. solution of the workload process of the system *)
        (* ====================================================== *)
        sM = Total[M[[k;;K]]];
		Qwmm = Sum[Dm[[i]],{i,k}];
        Qwpm = Table[0,{NN*sM}, {NN}];
        Qwmp = Table[0,{NN}, {NN*sM}];
        Qwpp = Table[0,{NN*sM},{NN*sM}];    
        kix = 1;
		Do[
            Qwmp[[All,kix;;kix+NN*M[[i]]-1]] = KroneckerProduct[Dm[[i+1]], {sigma[[i]]}];
            Qwpm[[kix;;kix+NN*M[[i]]-1,All]] = KroneckerProduct[II,s[[i]]];
            Qwpp[[kix;;kix+NN*M[[i]]-1,kix;;kix+NN*M[[i]]-1]] = KroneckerProduct[II,S[[i]]];
            kix += NN*M[[i]];
		,{i,k,K}];
		(* calculate fundamental matrices *)
        {Psiw, Kw, Uw} = FluidFundamentalMatrices[Qwpp, Qwpm, Qwmp, Qwmm, "PKU", prec];
        (* calculate boundary vector *)
        Ua = Table[1,{NN}] + 2*Total[Qwmp.Inverse[-Kw],{2}];   
		X = Join[Uw,Transpose[{Ua}],2];
		b = Table[0,{NN}];
		b = Append[b,1];
		{Qm,Rm}=QRDecomposition[Transpose[X]];
		pm=Flatten[Inverse[Rm].Qm.b];
		Bw = Table[0,{NN*sM},{NN}];
        Bw[[1;;NN*M[[k]],All]] = KroneckerProduct[II,s[[k]]];
        kappa = pm.Qwmp / Total[pm.Qwmp.Inverse[-Kw].Bw];
		
        If[k<K,
            (* step 2. construct fluid model for the remaining sojourn time process *)
            (* ==================================================================\[Equal] *)
            (* (for each class except the highest priority) *)
            Qsmm = Sum[Dm[[i]],{i,k+1}];
            Np = Dimensions[Kw][[1]];
            Qspm = Table[0,{Np+NN*Total[M[[k+1;;K]]]}, {NN}];
            Qsmp = Table[0,{NN}, {Np+NN*Total[M[[k+1;;K]]]}];
            Qspp = Table[0,{Np+NN*Total[M[[k+1;;K]]]},{Np+NN*Total[M[[k+1;;K]]]}];
            Qspp[[1;;Np,1;;Np]] = Kw;
            Qspm[[1;;Np,1;;NN]] = Bw;
            kix = Np+1;
			Do[
				Qsmp[[All,kix;;kix+NN*M[[i]]-1]] = KroneckerProduct[Dm[[i+1]], {sigma[[i]]}];
				Qspm[[kix;;kix+NN*M[[i]]-1,All]] = KroneckerProduct[II,s[[i]]];
				Qspp[[kix;;kix+NN*M[[i]]-1,kix;;kix+NN*M[[i]]-1]] = KroneckerProduct[II,S[[i]]];
				kix += NN*M[[i]];
			,{i,k+1,K}];         
			inis = Join[kappa, Table[0,{NN*Total[M[[k+1;;K]]]}]];
            Psis = FluidFundamentalMatrices[Qspp, Qspm, Qsmp, Qsmm, "P", prec];
			(* step 3. calculate the performance measures *)
            (* ========================================== *)
			argIx = 1;
			retIx = 1;
			While[argIx<=Length[varargin],
				res = {};
				If[MemberQ[eaten, argIx],
					argIx++;
					Continue[];
				,If[varargin[[argIx]]=="stMoms",
                    (* MOMENTS OF THE SOJOURN TIME *)
                    (* ~~~~~~~~~~~~~~~~~~~~~~~~~~~ *)
                    numOfSTMoms = varargin[[argIx+1]];
                    Pn = {Psis};
                    rtMoms = {};
					Do[
						A = Qspp + Psis.Qsmp;
                        B = Qsmm + Qsmp.Psis;
                        C = -2*n*Pn[[n]];
                        bino = 1;
						Do[
						    bino *= (n-i+1) / i;
                            C += bino * Pn[[i+1]].Qsmp.Pn[[n-i+1]];
                        ,{i,n-1}];
                        P = LyapunovSolve[A, B, -C];
                        AppendTo[Pn, P];
                        AppendTo[rtMoms,Total[inis.P*(-1)^n] / 2^n];
					,{n,numOfSTMoms}];
                    res = rtMoms;
                    argIx ++;
				,If[varargin[[argIx]]=="stDistr",
                    (* DISTRIBUTION OF THE SOJOURN TIME *)
                    (* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *)
                    stCdfPoints = varargin[[argIx+1]];
					Do[
					    L = erlMaxOrder;
                        lambda = L/t/2;
                        Psie = FluidFundamentalMatrices[Qspp-lambda*IdentityMatrix[Length[Qspp]], Qspm, Qsmp, Qsmm-lambda*IdentityMatrix[Length[Qsmm]], "P", prec];
                        Pn = {Psie};
                        pr = Total[inis.Psie];
						Do[
							A = Qspp + Psie.Qsmp - lambda*IdentityMatrix[Length[Qspp]];
                            B = Qsmm + Qsmp.Psie - lambda*IdentityMatrix[Length[Qsmm]];
                            C = 2*lambda*Pn[[n]];
							Do[C += Pn[[i+1]].Qsmp.Pn[[n-i+1]],{i,n-1}];
                            P = LyapunovSolve[A, B, -C];
                            AppendTo[Pn, P];
                            pr += Total[inis.P];
						,{n,L-1}];
                        AppendTo[res, pr];
                    ,{t,stCdfPoints}];
                    argIx++;
				,If[varargin[[argIx]]=="qlMoms",
                    (* MOMENTS OF THE NUMBER OF JOBS *)
                    (* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *)
					numOfQLMoms = varargin[[argIx+1]];
                    (* first calculate it at departure instants *)
                    QLDPn = {Psis};
                    dqlMoms = {};
                    Do[
                        A = Qspp + Psis.Qsmp;
                        B = Qsmm + Qsmp.Psis;
                        C = n*QLDPn[[n]].Dm[[k+1]];
                        bino = 1;
						Do[
                            bino *= (n-i+1) / i;
                            C += bino * QLDPn[[i+1]].Qsmp.QLDPn[[n-i+1]];
                        ,{i,n-1}];
                        P = LyapunovSolve[A, B, -C];
                        AppendTo[QLDPn, P];
                        AppendTo[dqlMoms, Total[inis.P]];
                    ,{n,numOfQLMoms}];
                    dqlMoms = MomsFromFactorialMoms[dqlMoms];                   
                    (* now calculate it at random time instance *)
                    pi = CTMCSolve[sD];
                    lambdak = Total[pi.Dm[[k+1]]];
                    QLPn = {pi};
                    iTerm = Inverse[Table[1,{NN},{1}].{pi} - sD];
					qlMoms = {};
					Do[
                        sumP = Total[inis.QLDPn[[n+1]]] + n*(inis.QLDPn[[n]] - QLPn[[n]].Dm[[k+1]]/lambdak).iTerm.Total[Dm[[k+1]],{2}];
                        P = sumP*pi + n*(QLPn[[n]].Dm[[k+1]] - inis.QLDPn[[n]]*lambdak).iTerm;
                        AppendTo[QLPn, P];
                        AppendTo[qlMoms, Total[P]];
                    ,{n,numOfQLMoms}];
                    qlMoms = MomsFromFactorialMoms[qlMoms];
                    res = qlMoms;
                    argIx++;
				,If[varargin[[argIx]]=="qlDistr",
                    (* DISTRIBUTION OF THE NUMBER OF JOBS *)
                    (* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *)
                    numOfQLProbs = varargin[[argIx+1]];
                    sDk = Sum[Dm[[i]],{i,k}];
                    (* first calculate it at departure instants *)
                    Psid = FluidFundamentalMatrices[Qspp, Qspm, Qsmp, sDk, "P", prec];
                    Pn = {Psid};
                    dqlProbs = {inis.Psid};
					Do[
                        A = Qspp + Psid.Qsmp;
                        B = sDk + Qsmp.Psid;
                        C = Pn[[n]].Dm[[k+1]];
						Do[C+=Pn[[i+1]].Qsmp.Pn[[n-i+1]],{i,n-1}];
                        P = LyapunovSolve[A, B, -C];
                        AppendTo[Pn, P];
                        AppendTo[dqlProbs, inis.P];
                    ,{n,numOfQLProbs-1}];
                    (* now calculate it at random time instance *)
                    pi = CTMCSolve[sD];
                    lambdak = Total[pi.Dm[[k+1]]];
                    iTerm = Inverse[-(sD-Dm[[k+1]])];
                    qlProbs = {lambdak*dqlProbs[[1,All]].iTerm};
                    Do[
                        P = (qlProbs[[n,All]].Dm[[k+1]]+lambdak*(dqlProbs[[n+1,All]]-dqlProbs[[n,All]])).iTerm;
                        AppendTo[qlProbs, P];
                    ,{n,numOfQLProbs-1}];
                    qlProbs = Total[qlProbs,{2}];
                    res = qlProbs;		
                    argIx++;
		        ,
					Throw[StringJoin["MMAPPH1PRPR: Unknown parameter ", ToString[varargin[[argIx]]]]];
				]]]]];
				If[retIx>Length[Ret],AppendTo[Ret,Transpose[{res}]],Ret[[retIx]]=Join[Ret[[retIx]],Transpose[{res}],2]];
				retIx++;
				argIx++;
			];
		,
            (* step 3. calculate the performance measures *)
            (* ========================================== *)
            retIx = 1;
            argIx = 1;
			While[argIx<=Length[varargin],
				res = {};
				If[MemberQ[eaten, argIx],
					argIx++;
					Continue[];
				,If[varargin[[argIx]]=="stMoms",
                    (* MOMENTS OF THE SOJOURN TIME *)
                    (* ~~~~~~~~~~~~~~~~~~~~~~~~~~~ *)
                    numOfSTMoms = varargin[[argIx+1]];
                    res = Table[i!*kappa.MatrixPower[Inverse[-Kw],i+1].Total[Bw,{2}],{i,numOfSTMoms}];
                    argIx++;
				,If[varargin[[argIx]]=="stDistr",
                    (* DISTRIBUTION OF THE SOJOURN TIME *)
                    (* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *)
					stCdfPoints = varargin[[argIx+1]];
                    res=Table[kappa.Inverse[-Kw].(IdentityMatrix[Length[Kw]]-MatrixExp[Kw*t]).Total[Bw,{2}],{t,stCdfPoints}];
                    argIx++;
				,If[varargin[[argIx]]=="qlMoms" || varargin[[argIx]]=="qlDistr",
			        L = KroneckerProduct[sD-Dm[[k+1]],IdentityMatrix[M[[k]]]]+KroneckerProduct[IdentityMatrix[NN],S[[k]]];
                    B = KroneckerProduct[IdentityMatrix[NN],s[[k]].{sigma[[k]]}];
                    F = KroneckerProduct[Dm[[k+1]],IdentityMatrix[M[[k]]]];
                    L0 = KroneckerProduct[sD-Dm[[k+1]],IdentityMatrix[M[[k]]]];
                    R = QBDFundamentalMatrices[B, L, F, "R", prec];
                    p0 = CTMCSolve[L0+R.B];
                    p0 = p0/Total[p0.Inverse[IdentityMatrix[Length[R]]-R]];
					If[varargin[[argIx]]=="qlMoms",
                        (* MOMENTS OF THE NUMBER OF JOBS *)
                        (* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *)
                        numOfQLMoms = varargin[[argIx+1]];
                        res = MomsFromFactorialMoms[Table[Total[i!*p0.MatrixPower[R,i].MatrixPower[Inverse[IdentityMatrix[Length[R]]-R],i+1]],{i,numOfQLMoms}]];
					, If[varargin[[argIx]]=="qlDistr",
                        qlProbs = {p0};
						Rpow = R;
                        Do[AppendTo[qlProbs,p0.Rpow]; Rpow=Rpow.R,{i,numOfQLProbs-1}];
                        res = Total[qlProbs,{2}];
					]];
					argIx++;
		        ,
					Throw[StringJoin["MMAPPH1PRPR: Unknown parameter ", ToString[varargin[[argIx]]]]];
				]]]];
				If[retIx>Length[Ret],AppendTo[Ret,Transpose[{res}]],Ret[[retIx]]=Join[Ret[[retIx]],Transpose[{res}],2]];
				retIx++;
				argIx++;
			];
		];
	,{k,classes}];
	Return[Ret];
];


MMAPPH1NPPR[Dm_,sigma_,S_,argv__]:=
Module[{varargin,prec,eaten,classes,K,erlMaxOrder,
D0,s,sD,M,NN,II,Ret,sM,Qwmm,Qwmp,Qwpm,Qwpp,Qwzp,
Qwmz,Qwpz,Qwzz,Qwzm,Mlo,Mhi,ro,bs,kix,Dlo,kix2,bs2,Psikw,
Qkwpp,Qkwpm,Qkwpz,Qkwzp,Qkwzz,Qkwzm,Qkwmp,Qkwmz,Qkwmm,
Psiw,Kw,Uw,Ua,b,X,pm,Qm,Rm,Bw,kappa,argIx,
numOfSTMoms,res,A,B,C,bino,Qsmm,Qsmp,Qspm,Qspp,Np,
inis,Psis,P,Pn,rtMoms,stCdfPoints,pr,lambda,L,Psie,
numOfQLMoms,QLDPn,dqlMoms,lambdak,pi,QLPn,iTerm,sumP,
qlMoms,numOfQLProbs,sDk,Psid,dqlProbs,qlProbs,retIx,
F,L0,p0,R,Rpow,lambdaS,phi,q0,qL,pk,Qwzpk,vix,Bk,ztag,
Mx,qLkp1,qLii,Ak,V1,Qwmpk,sD0k,BM,CM,DM,Kwu,Bwu,iniw,pwu,norm,
AM,KN,wtMoms,Pnr,lambdae,W,iW,w,omega,Psii,XDn,z,Z,zeta,Tmp},
	varargin = List[argv];
	K = Length[Dm]-1;
    (* parse options *)
    prec = N[10^-14];
    erlMaxOrder = 200;
	classes = Range[K];
    eaten = {};
	Do[
		If[varargin[[i]]=="prec",
			prec = varargin[[i+1]];
			AppendTo[eaten,i]; 
			AppendTo[eaten,i+1];
		, If[varargin[[i]]=="erlMaxOrder",
			erlMaxOrder = varargin[[i+1]];
			AppendTo[eaten,i]; 
			AppendTo[eaten,i+1];
		, If[varargin[[i]]=="classes",
			classes = varargin[[i+1]];
			AppendTo[eaten,i]; 
			AppendTo[eaten,i+1];
		]]];
	,{i,Length[varargin]}];

    If[BuTools`CheckInput && Not[CheckMMAPRepresentation[Dm]], Throw["MMAPPH1NPPR: The arrival process is not a valid MMAP representation!"]]; 
    Do[If[BuTools`CheckInput && Not[CheckPHRepresentation[sigma[[k]],S[[k]]]], Throw["MMAPPH1NPPR: the vector and matrix describing the service times is not a valid PH representation!"]],{k,K}];

    (* some preparation *)
    D0 = Dm[[1]];
    NN = Dimensions[D0][[1]];
    II = IdentityMatrix[NN];
    sD = Sum[Dm[[i]],{i,K+1}];
	s = Table[Transpose[{Total[-S[[i]],{2}]}],{i,K}];
	M = Table[Length[sigma[[i]]],{i,K}];
	
    (* step 1. solution of the workload process of the joint queue *)
    (* =========================================================\[Equal] *)
    Qwmm = D0;
    Qwpp = Table[0,{NN*Total[M]},{NN*Total[M]}];
    Qwmp = Table[0,{NN},{NN*Total[M]}];
    Qwpm = Table[0,{NN*Total[M]},{NN}];
    kix = 1;
	Do[
		bs = NN*M[[i]];
        Qwpp[[kix;;kix+bs-1,kix;;kix+bs-1]] = KroneckerProduct[II,S[[i]]];
        Qwmp[[All,kix;;kix+bs-1]] = KroneckerProduct[Dm[[i+1]], {sigma[[i]]}];
        Qwpm[[kix;;kix+bs-1,All]] = KroneckerProduct[II,s[[i]]];
        kix += bs;
	,{i,K}];
	(* calculate fundamental matrices *)
    {Psiw, Kw, Uw} = FluidFundamentalMatrices[Qwpp, Qwpm, Qwmp, Qwmm, "PKU", prec];
	(* calculate boundary vector *)
    Ua = Table[1,{NN}] + 2*Total[Qwmp.Inverse[-Kw],{2}];   
	X = Join[Uw,Transpose[{Ua}],2];
	b = Table[0,{NN}];
	b = Append[b,1];
	{Qm,Rm}=QRDecomposition[Transpose[X]];
	pm=Flatten[Inverse[Rm].Qm.b];
    ro = ((1-Total[pm])/2)/(Total[pm]+(1-Total[pm])/2); (* calc idle time with weight=1, and the busy time with weight=1/2 *)
    kappa = pm/Total[pm];
      
    pi = CTMCSolve[sD];
	lambda = Table[Total[pi.Dm[[i+1]]],{i,K}];
    
	Psiw={};
	Qwmp = {}; Qwzp = {}; Qwpp = {};
	Qwmz = {}; Qwpz = {}; Qwzz = {};
	Qwmm = {}; Qwpm = {}; Qwzm = {};
    Do[
        (* step 2. construct a workload process for classes k...K *)
        (* ====================================================== *)
        Mlo = Total[M[[1;;k-1]]];
        Mhi = Total[M[[k;;K]]];

        Qkwpp = Table[0,{NN*Mlo*Mhi+NN*Mhi},{NN*Mlo*Mhi+NN*Mhi}];
        Qkwpz = Table[0,{NN*Mlo*Mhi+NN*Mhi},{NN*Mlo}]; 
        Qkwpm = Table[0,{NN*Mlo*Mhi+NN*Mhi}, {NN}];
        Qkwmz = Table[0,{NN}, {NN*Mlo}];
        Qkwmp = Table[0,{NN}, {NN*Mlo*Mhi+NN*Mhi}];
        Dlo = Sum[Dm[[i]],{i,k}];
        Qkwmm = Dlo;
        Qkwzp = Table[0,{NN*Mlo}, {NN*Mlo*Mhi+NN*Mhi}];
        Qkwzm = Table[0,{NN*Mlo},{NN}];
        Qkwzz = Table[0,{NN*Mlo},{NN*Mlo}];
        kix = 1;
		Do[
            kix2 = 1;
            Do[
                bs = NN*M[[j]]*M[[i]];
                bs2 = NN*M[[j]];
                Qkwpp[[kix;;kix+bs-1,kix;;kix+bs-1]] = KroneckerProduct[II,KroneckerProduct[IdentityMatrix[M[[j]]],S[[i]]]];
                Qkwpz[[kix;;kix+bs-1,kix2;;kix2+bs2-1]] = KroneckerProduct[II,KroneckerProduct[IdentityMatrix[M[[j]]],s[[i]]]];
                Qkwzp[[kix2;;kix2+bs2-1,kix;;kix+bs-1]] = KroneckerProduct[Dm[[i+1]],KroneckerProduct[IdentityMatrix[M[[j]]], {sigma[[i]]}]];
                kix += bs;
                kix2 += bs2;
            ,{j,k-1}]
        ,{i,k,K}];
        Do[
            bs = NN*M[[i]];
            Qkwpp[[kix;;kix+bs-1,kix;;kix+bs-1]] = KroneckerProduct[II,S[[i]]];
            Qkwpm[[kix;;kix+bs-1,All]] = KroneckerProduct[II,s[[i]]];
            Qkwmp[[All,kix;;kix+bs-1]] = KroneckerProduct[Dm[[i+1]],{sigma[[i]]}];
            kix += bs;
        ,{i,k,K}];
        kix = 1;
        Do[
            bs = NN*M[[j]];
            Qkwzz[[kix;;kix+bs-1,kix;;kix+bs-1]] = KroneckerProduct[Dlo, IdentityMatrix[M[[j]]]] + KroneckerProduct[II, S[[j]]];
            Qkwzm[[kix;;kix+bs-1,All]] = KroneckerProduct[II, s[[j]]];
            kix += bs;
        ,{j,k-1}];
        Psikw = FluidFundamentalMatrices[If[Mlo==0,Qkwpp,Qkwpp+Qkwpz.Inverse[-Qkwzz].Qkwzp], If[Mlo==0,Qkwpm,Qkwpm+Qkwpz.Inverse[-Qkwzz].Qkwzm], Qkwmp, Qkwmm, "P", prec];
        AppendTo[Psiw, Psikw];
        AppendTo[Qwzp, Qkwzp]; AppendTo[Qwmp, Qkwmp]; AppendTo[Qwpp, Qkwpp]; 
        AppendTo[Qwmz, Qkwmz]; AppendTo[Qwpz, Qkwpz]; AppendTo[Qwzz, Qkwzz]; 
		AppendTo[Qwmm, Qkwmm]; AppendTo[Qwpm, Qkwpm]; AppendTo[Qwzm, Qkwzm]; 
    ,{k,K}];
    (* step 3. calculate Phi vectors *)
    (* ===========================\[Equal] *)
    lambdaS = Total[lambda];
    phi = {(1-ro)*kappa.(-D0) / lambdaS};
    q0 = {{}};
    qL = {{}};
    Do[
		sDk = Sum[Dm[[j]],{j,k+1}];
        (* pk *)
        pk = Total[lambda[[1;;k]]]/lambdaS - (1-ro)*kappa.Total[sDk,{2}]/lambdaS;
        (* A^(k,1) *)
        Qwzpk = Qwzp[[k+1]];
        vix = 1;
		Ak = {};
        Do[
            bs = NN*M[[ii]];
            V1 = Qwzpk[[vix;;vix+bs-1,All]];
            AppendTo[Ak, KroneckerProduct[II,{sigma[[ii]]}].Inverse[-KroneckerProduct[sDk,IdentityMatrix[M[[ii]]]]-KroneckerProduct[II,S[[ii]]]].(KroneckerProduct[II,s[[ii]]] + V1.Psiw[[k+1]])];
            vix += bs;
        ,{ii,k}];
        (* B^k *)
        Qwmpk = Qwmp[[k+1]];
        Bk = Qwmpk.Psiw[[k+1]];
        ztag = phi[[1]].(Inverse[-D0].Dm[[k+1]].Ak[[k]] - Ak[[1]] + Inverse[-D0].Bk);
        Do[
            ztag += phi[[i+1]].(Ak[[i]]-Ak[[i+1]]) + phi[[1]].Inverse[-D0].Dm[[i+1]].Ak[[i]];
        ,{i,k-1}];
        Mx = IdentityMatrix[Length[Ak[[k]]]]-Ak[[k]];
        Mx[[All,1;;1]] = Table[1,{NN},{1}];
        AppendTo[phi, Prepend[ztag[[2;;]],pk].Inverse[Mx]]; (* phi(k) = Psi^(k)_k * p(k). Psi^(k)_i = phi(i) / p(k)] *)
        AppendTo[q0, phi[[1]].Inverse[-D0]];
        qLkp1 = {};
        Do[
            qLii = (phi[[ii+1]] - phi[[ii]] + phi[[1]].Inverse[-D0].Dm[[ii+1]]).KroneckerProduct[II,{sigma[[ii]]}].Inverse[-KroneckerProduct[sDk,IdentityMatrix[M[[ii]]]]-KroneckerProduct[II,S[[ii]]]];
            qLkp1 = Join[qLkp1, qLii];
        ,{ii,k}];
		AppendTo[qL, qLkp1];
    ,{k,K-1}];

    (* step 4. calculate performance measures *)
    (* ====================================\[Equal] *)
    Ret = {};
	Do[
		sD0k = Sum[Dm[[i]],{i,k}];
		If[k<K,
            (* step 4.1 calculate distribution of the workload process right  *)
            (* before the arrivals of class k jobs *)
            (* ============================================================ *)
            Kw = If[Length[Qwzz[[k]]]==0,Qwpp[[k]],Qwpp[[k]]+Qwpz[[k]].Inverse[-Qwzz[[k]]].Qwzp[[k]]] + Psiw[[k]].Qwmp[[k]];                   
            BM = {{}}; CM = {}; DM = {{}};
            Do[
				Tmp = KroneckerProduct[II,S[[i]]];
                BM = If[i==1,Tmp,ArrayFlatten[{{BM,0},{0,Tmp}}]];
                CM = Join[CM, KroneckerProduct[II,s[[i]]]];
				Tmp = KroneckerProduct[Dm[[k+1]],IdentityMatrix[M[[i]]]];
                DM = If[i==1,Tmp,ArrayFlatten[{{DM,0},{0,Tmp}}]];
            ,{i,k-1}];
            Kwu = If[Length[Qwzz[[k]]]==0,Kw,ArrayFlatten[{{Kw,(Qwpz[[k]]+Psiw[[k]].Qwmz[[k]]).Inverse[-Qwzz[[k]]].DM},{Table[0,{Dimensions[BM][[1]]},{Dimensions[Kw][[2]]}], BM}}]];
            Bwu = Join[Psiw[[k]].Dm[[k+1]], CM];
            If[k>1,
			    iniw = Join[q0[[k]].Qwmp[[k]]+qL[[k]].Qwzp[[k]], qL[[k]].DM];
                pwu = q0[[k]].Dm[[k+1]];
            ,
                iniw = pm.Qwmp[[k]];
                pwu = pm.Dm[[k+1]];
            ];
            norm = Total[pwu] + Total[iniw.Inverse[-Kwu].Bwu];
            pwu /= norm;
            iniw /= norm;
			(* step 4.2 create the fluid model whose first passage time equals the *)
            (* WAITING time of the low prioroity customers *)
            (* ================================================================== *)
            KN = Dimensions[Kwu][[1]];
            Qspp = Table[0,{KN+NN*Total[M[[k+1;;]]]},{KN+NN*Total[M[[k+1;;]]]}];
            Qspm = Table[0,{KN+NN*Total[M[[k+1;;]]]}, {NN}];
            Qsmp = Table[0,{NN}, {KN+NN*Total[M[[k+1;;]]]}];
            Qsmm = sD0k + Dm[[k+1]];
            kix = 1;
            Do[
                bs = NN*M[[i]];
                Qspp[[KN+kix;;KN+kix+bs-1,KN+kix;;KN+kix+bs-1]] = KroneckerProduct[II,S[[i]]];
                Qspm[[KN+kix;;KN+kix+bs-1,All]] = KroneckerProduct[II,s[[i]]];
                Qsmp[[All,KN+kix;;KN+kix+bs-1]] = KroneckerProduct[Dm[[i+1]],{sigma[[i]]}];
                kix += bs;
            ,{i,k+1,K}];
            Qspp[[1;;KN,1;;KN]] = Kwu;
            Qspm[[1;;KN,All]] = Bwu;           
            inis = Join[iniw, Table[0,{NN*Total[M[[k+1;;]]]}]];
            (* calculate fundamental matrix *)
            Psis = FluidFundamentalMatrices[Qspp, Qspm, Qsmp, Qsmm, "P", prec];
            
			(* step 4.3. calculate the performance measures *)
            (* ========================================== *)
			argIx = 1;
			retIx = 1;
			While[argIx<=Length[varargin],
				res = {};
				If[MemberQ[eaten, argIx],
					argIx++;
					Continue[];
				,If[varargin[[argIx]]=="stMoms",
                    (* MOMENTS OF THE SOJOURN TIME *)
                    (* ~~~~~~~~~~~~~~~~~~~~~~~~~~~ *)
                    numOfSTMoms = varargin[[argIx+1]];
                    Pn = {Psis};
                    wtMoms = {};
					Do[
						A = Qspp + Psis.Qsmp;
                        B = Qsmm + Qsmp.Psis;
                        C = -2*n*Pn[[n]];
                        bino = 1;
						Do[
						    bino *= (n-i+1) / i;
                            C += bino * Pn[[i+1]].Qsmp.Pn[[n-i+1]];
                        ,{i,n-1}];
                        P = LyapunovSolve[A, B, -C];
                        AppendTo[Pn, P];
                        AppendTo[wtMoms,Total[inis.P*(-1)^n] / 2^n];
					,{n,numOfSTMoms}];
                    (* calculate RESPONSE time moments *)
                    Pnr = {Total[inis.Pn[[1]]]*sigma[[k]]};
                    rtMoms = {};
                    Do[
                        P = n*Pnr[[n]].Inverse[-S[[k]]] + (-1)^n*Total[inis.Pn[[n+1]]]*sigma[[k]] / 2^n;
                        AppendTo[Pnr,P];
                        AppendTo[rtMoms, Total[P]+Total[pwu]*n!*Total[sigma[[k]].MatrixPower[Inverse[-S[[k]]],n]]];
                    ,{n,numOfSTMoms}];
                    res = rtMoms;
                    argIx ++;
				,If[varargin[[argIx]]=="stDistr",
                    (* DISTRIBUTION OF THE SOJOURN TIME *)
                    (* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *)
                    stCdfPoints = varargin[[argIx+1]];
					Do[
					    L = erlMaxOrder;
                        lambdae = L/t/2;
                        Psie = FluidFundamentalMatrices[Qspp-lambdae*IdentityMatrix[Length[Qspp]], Qspm, Qsmp, Qsmm-lambdae*IdentityMatrix[Length[Qsmm]], "P", prec];
                        Pn = {Psie};
						pr = (Total[pwu] + Total[inis.Psie]) * (1-Total[sigma[[k]].MatrixPower[Inverse[IdentityMatrix[Length[S[[k]]]]-S[[k]]/2/lambdae],L]]);
                        Do[
							A = Qspp + Psie.Qsmp - lambdae*IdentityMatrix[Length[Qspp]];
                            B = Qsmm + Qsmp.Psie - lambdae*IdentityMatrix[Length[Qsmm]];
                            C = 2*lambdae*Pn[[n]];
							Do[C += Pn[[i+1]].Qsmp.Pn[[n-i+1]],{i,n-1}];
                            P = LyapunovSolve[A, B, -C];
                            AppendTo[Pn, P];
                            pr += Total[inis.P] * (1-Total[sigma[[k]].MatrixPower[Inverse[IdentityMatrix[Length[S[[k]]]]-S[[k]]/2/lambdae],L-n]]);
						,{n,L-1}];
                        AppendTo[res, pr];
                    ,{t,stCdfPoints}];
                    argIx++;
				,If[varargin[[argIx]]=="qlMoms" || varargin[[argIx]]=="qlDistr",
                    W = Inverse[-KroneckerProduct[sD-Dm[[k+1]],IdentityMatrix[M[[k]]]]-KroneckerProduct[II,S[[k]]]].KroneckerProduct[Dm[[k+1]],IdentityMatrix[M[[k]]]];
                    iW = Inverse[IdentityMatrix[Length[W]]-W];
                    w = KroneckerProduct[II,{sigma[[k]]}];
                    omega = Inverse[-KroneckerProduct[sD-Dm[[k+1]],IdentityMatrix[M[[k]]]]-KroneckerProduct[II,S[[k]]]].KroneckerProduct[II,s[[k]]];
					If[varargin[[argIx]]=="qlMoms",
						(* MOMENTS OF THE NUMBER OF JOBS *)
						(* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *)
						numOfQLMoms = varargin[[argIx+1]];
						(* first calculate it at departure instants *)
						Psii = {Psis};
						QLDPn = {inis.Psii[[1]].w.iW};
						Do[
							A = Qspp + Psis.Qsmp;
							B = Qsmm + Qsmp.Psis;
							C = n*Psii[[n]].Dm[[k+1]];
							bino = 1;
							Do[
								bino *= (n-i+1) / i;
								C += bino * Psii[[i+1]].Qsmp.Psii[[n-i+1]];
							,{i,n-1}];
							P = LyapunovSolve[A, B, -C];
							AppendTo[Psii, P];
							AppendTo[QLDPn, n*QLDPn[[n]].iW.W + inis.P.w.iW];
						,{n,numOfQLMoms}];
						Do[
							QLDPn[[n+1]] = (QLDPn[[n+1]]+pwu.w.MatrixPower[iW,n+1].If[n==0,IdentityMatrix[Length[W]],MatrixPower[W,n]]).omega;
                        ,{n,0,numOfQLMoms}];
						(* now calculate it at random time instance *)
						QLPn = {pi};
						iTerm = Inverse[Table[1,{NN},{1}].{pi} - sD];
						qlMoms = {};
						Do[
							sumP = Total[QLDPn[[n+1]]] + n*(QLDPn[[n]] - QLPn[[n]].Dm[[k+1]]/lambda[[k]]).iTerm.Total[Dm[[k+1]],{2}];
							P = sumP*pi + n*(QLPn[[n]].Dm[[k+1]] - QLDPn[[n]]*lambda[[k]]).iTerm;
							AppendTo[QLPn, P];
							AppendTo[qlMoms, Total[P]];
						,{n,numOfQLMoms}];
						qlMoms = MomsFromFactorialMoms[qlMoms];
						res = qlMoms;
						argIx++;
					,If[varargin[[argIx]]=="qlDistr",
						(* DISTRIBUTION OF THE NUMBER OF JOBS *)
						(* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *)
						numOfQLProbs = varargin[[argIx+1]];
						(* first calculate it at departure instants *)
						Psid = FluidFundamentalMatrices[Qspp, Qspm, Qsmp, sD0k, "P", prec];
						Pn = {Psid};
						XDn = inis.Psid.w;
						dqlProbs = {(XDn+pwu.w).omega};
						Do[
							A = Qspp + Psid.Qsmp;
							B = sD0k + Qsmp.Psid;
							C = Pn[[n]].Dm[[k+1]];
							Do[C+=Pn[[i+1]].Qsmp.Pn[[n-i+1]],{i,n-1}];
							P = LyapunovSolve[A, B, -C];
							AppendTo[Pn, P];
							XDn = XDn.W + inis.P.w;
							AppendTo[dqlProbs, (XDn+pwu.w.MatrixPower[W,n]).omega];
						,{n,numOfQLProbs-1}];
						(* now calculate it at random time instance *)
						iTerm = Inverse[-(sD-Dm[[k+1]])];
						qlProbs = {lambda[[k]]*dqlProbs[[1,All]].iTerm};
						Do[
							P = (qlProbs[[n,All]].Dm[[k+1]]+lambda[[k]]*(dqlProbs[[n+1,All]]-dqlProbs[[n,All]])).iTerm;
							AppendTo[qlProbs, P];
						,{n,numOfQLProbs-1}];
						qlProbs = Total[qlProbs,{2}];
						res = qlProbs;		
						argIx++;
					]];
		        ,
					Throw[StringJoin["MMAPPH1NPPR: Unknown parameter ", ToString[varargin[[argIx]]]]];
				]]]];
				If[retIx>Length[Ret],AppendTo[Ret,Transpose[{res}]],Ret[[retIx]]=Join[Ret[[retIx]],Transpose[{res}],2]];
				retIx++;
				argIx++;
			];
		,
            (* step 3. calculate the performance measures *)
            (* ========================================== *)
            retIx = 1;
            argIx = 1;
			While[argIx<=Length[varargin],
				res = {};
				If[MemberQ[eaten, argIx],
					argIx++;
					Continue[];
				, If[varargin[[argIx]]=="stMoms" || varargin[[argIx]]=="stDistr",
                    Kw = Qwpp[[k]]+Qwpz[[k]].Inverse[-Qwzz[[k]]].Qwzp[[k]] + Psiw[[k]].Qwmp[[k]];                   
                    AM = {{}}; BM = {{}}; CM = {}; DM = {{}};
                    Do[
                        Tmp = KroneckerProduct[Table[1,{NN},{1}],KroneckerProduct[IdentityMatrix[M[[i]]],s[[k]]]];
						AM = If[i==1,Tmp,ArrayFlatten[{{AM,0},{0,Tmp}}]];
                        BM = If[i==1,S[[i]],ArrayFlatten[{{BM,0},{0,S[[i]]}}]];
                        CM = Join[CM, s[[i]]];
						Tmp = KroneckerProduct[Dm[[k+1]],IdentityMatrix[M[[i]]]];
                        DM = If[i==1,Tmp,ArrayFlatten[{{DM,0},{0,Tmp}}]];
                    ,{i,k-1}];
                    Z = ArrayFlatten[{{Kw, Join[AM,Table[0,{NN*M[[k]]},{Dimensions[AM][[2]]}]]},{Table[0,{Dimensions[BM][[1]]},{Dimensions[Kw][[2]]}], BM}}];
                    z = ArrayFlatten[{{Table[0,{Dimensions[AM][[1]]},{1}]},{KroneckerProduct[Table[1,{NN},{1}],s[[k]]]},{CM}}];
                    iniw = Join[q0[[k]].Qwmp[[k]]+qL[[k]].Qwzp[[k]], Table[0,{Length[BM]}]];
                    zeta = iniw/Total[iniw.Inverse[-Z].z];   
					If[varargin[[argIx]]=="stMoms",
					    numOfSTMoms = varargin[[argIx+1]];
						res = Table[i!*(zeta.MatrixPower[Inverse[-Z],(i+1)].z)[[1]],{i,numOfSTMoms}];
                    ,
						stCdfPoints = varargin[[argIx+1]];
						res = Table[(zeta.Inverse[-Z].(IdentityMatrix[Length[Z]]-MatrixExp[Z*t]).z)[[1]],{t,stCdfPoints}];
                    ];
					argIx++;
				, If[varargin[[argIx]]=="qlMoms" || varargin[[argIx]]=="qlDistr",
                    L = Table[0,{NN*Total[M]},{NN*Total[M]}];
                    B = Table[0,{NN*Total[M]},{NN*Total[M]}];
                    F = Table[0,{NN*Total[M]},{NN*Total[M]}];
                    kix = 1;
                    Do[
                        bs = NN*M[[i]];
                        F[[kix;;kix+bs-1,kix;;kix+bs-1]] = KroneckerProduct[Dm[[k+1]],IdentityMatrix[M[[i]]]];
                        L[[kix;;kix+bs-1,kix;;kix+bs-1]] = KroneckerProduct[sD0k,IdentityMatrix[M[[i]]]] + KroneckerProduct[II,S[[i]]];
                        If[i<K,
                            L[[kix;;kix+bs-1,NN*Total[M[[1;;k-1]]]+1;;]] = KroneckerProduct[II,s[[i]].{sigma[[k]]}];
                        ,
                            B[[kix;;kix+bs-1,NN*Total[M[[1;;k-1]]]+1;;]] = KroneckerProduct[II,s[[i]].{sigma[[k]]}];
                        ];
                        kix += bs;
                    ,{i,K}];
			        R = QBDFundamentalMatrices[B, L, F, "R", prec];
                    p0 = Join[qL[[k]], q0[[k]].KroneckerProduct[II,{sigma[[k]]}]];
                    p0 = p0/Total[p0.Inverse[IdentityMatrix[Length[R]]-R]];
					If[varargin[[argIx]]=="qlMoms",
						numOfQLMoms = varargin[[argIx+1]];
						res = MomsFromFactorialMoms[Table[Total[i!*p0.MatrixPower[R,i].MatrixPower[Inverse[IdentityMatrix[Length[R]]-R],i+1]],{i,numOfQLMoms}]];
					,
						numOfQLProbs = varargin[[argIx+1]];
                        qlProbs = Table[p0.If[i==0,IdentityMatrix[Length[R]],MatrixPower[R,i]],{i,0,numOfQLProbs-1}];
                        res = Total[qlProbs,{2}];  
					];
					argIx++;
		        ,
					Throw[StringJoin["MMAPPH1NPPR: Unknown parameter ", ToString[varargin[[argIx]]]]];
				]]];
				If[retIx>Length[Ret],AppendTo[Ret,Transpose[{res}]],Ret[[retIx]]=Join[Ret[[retIx]],Transpose[{res}],2]];
				retIx++;
				argIx++;
			];
		];
	,{k,classes}];
	Return[Ret];
];


End[(* Private *)];
EndPackage[];
