(* ::Package:: *)

(*
   BuTools Queues Package
*)

BeginPackage["BuTools`Queues`"];
QBDQueue::usage="";
MAPMAP1::usage="";
FluidQueue::usage="";
FluFluQueue::usage="";


Begin["`Private`"];


Needs["BuTools`MC`"];
Needs["BuTools`MAM`"];
Needs["BuTools`RepTrans`"];
Needs["BuTools`Moments`"];
Needs["BuTools`MAP`"];
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
		iniKi=Flatten[Inverse[Rm].Qm.-ini];
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
			iniKi=Flatten[Inverse[Rm].Qm.-ini];
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
			iniKi=Flatten[Inverse[Rm].Qm.-ini];
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
			iniKi=Flatten[Inverse[Rm].Qm.-inih];
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


End[(* Private *)];
EndPackage[];
