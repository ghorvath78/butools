(* ::Package:: *)

(*
   BuTools MAM Package
*)

BeginPackage["BuTools`MAM`"];
QBDFundamentalMatrices::usage = "M = QBDFundamentalMatrices[B, L, F, matrices, precision, maxNumIt, method]: Returns any combination of the R, G, and U matrices of a discrete or continuous time QBD";
QBDSolve::usage = "{pi0, R} = QBDSolve [B, L, F, L0, prec]: Returns the parameters of the matrix-geometrically distributed stationary distribution of a QBD";
QBDStationaryDistr::usage = "pi = QBDStationaryDistr [pi0, R, K]: Returns the stationary distribution of a QBD up to a given level";
FluidFundamentalMatrices::usage = "M = FluidFundamentalMatrices[Fpp, Fpm, Fmp, Fmm, matrices, precision, maxNumIt, method]: Returns any combination of the Psi, K, and U matrices of a canonical Markovian fluid model";
FluidSolve::usage = "{mass0, ini, K, clo} = FluidSolve [Fpp, Fpm, Fmp, Fmm, prec]: Returns the parameters of the matrix-exponentially distributed stationary distribution of a canonical Markovian fluid model";
FluidStationaryDistr::usage = "pi = FluidStationaryDistr [mass0, ini, K, clo, x]: Returns the cummulative distribution of the stationary fluid level of a Markovian fluid model";
GeneralFluidSolve::usage = "{mass0, ini, K, clo} = GeneralFluidSolve [Q, R, Q0, prec]: Returns the parameters of the matrix-exponentially distributed stationary distribution of a general Markovian fluid model";
MG1FundamentalMatrix::usage = "G = MG1FundamentalMatrix[A, precision, maxNumIt, method]: Returns matrix G, the fundamental matrix of an M/G/1 type Markov chain";
MG1StationaryDistr::usage = "pi = MG1StationaryDistr [A, B, G, K]: Returns the stationary distribution of an M/G/1 type Markov chain up to a given level";
GM1FundamentalMatrix::usage = "R = GM1FundamentalMatrix[A, precision, maxNumIt, method]: Returns matrix R, the fundamental matrix of a G/M/1 type Markov chain";
GM1StationaryDistr::usage = "pi = GM1StationaryDistr [B, R, K]: Returns the stationary distribution of a G/M/1 type Markov chain up to a given level";


Begin["`Private`"];


Needs["BuTools`MC`"];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckInput"]],BuTools`CheckInput=True];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckPrecision"]],BuTools`CheckPrecision=N[10^-12]];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`Verbose"]],BuTools`Verbose=False];


QBDFundamentalMatrices[B_, L_, F_, matrices_:"G", precision_:N[10^-14], maxNumIt_:50, method_:"CR", shift_:True]:=
Module[{m, II,continuous,lamb,Bm,Lm,Fm,theta,drift,Lold,uT,Bold,Fold,
BF,BB,G,PI,check,numit,Lstar,Bstar,Fstar,ret,U,R,resNorm},
    m = Length[L];
    II = IdentityMatrix[m];
    If[method!="CR" && BuTools`Verbose,Print["Warning: Currently only the 'CR' method is available in the Python implementation!"]];
    (* Convert to discrete time problem, if needed *)
    continuous = False;
    If[AnyTrue[Diagonal[L],Negative],
        continuous = True;
        lamb = Max[-Diagonal[L]];
        Bm = B / lamb; Lm = L / lamb + II; Fm = F / lamb;
	,
		Bm=B; Lm=L; Fm=F];
    
    (* Shift technique *)
    If[shift,
        theta = DTMCSolve[Bm+Lm+Fm];
        drift = Total[theta.Bm] - Total[theta.Fm];
        Lold = Lm;
        If[drift < 0, (* MC is transient -> use the dual MC*)
            Fold = Fm;
            Fm = Fm - ConstantArray[1,{m,1}].{theta.Fm};
            Lm = Lm + ConstantArray[1,{m,1}].{theta.Bm};
        ,
            uT = ConstantArray[1,{1,m}] / m;
            Bold = Bm;
            Bm = Bm - Transpose[{Total[Bm,{2}]}].uT;
            Lm = Lm + Transpose[{Total[Fm,{2}]}].uT;
		];
	];
    (* Start of Logaritmic Reduction (Basic) *)
    BF = Inverse[II - Lm];
    BB = BF.Fm;
    BF = BF.Bm;
    G = BF;
    PI = BB;
    check = 1;
    numit = 0;
    While[check > precision && numit < maxNumIt,
        Lstar = BF.BB + BB.BF;
        Bstar = BB.BB;
        Fstar = BF.BF;
        BB = Inverse[II-Lstar];
        BF = BB.Fstar;
        BB = BB.Bstar;
        G += PI.BF;
        PI = PI.BB;
        check = Min[Norm[BB,Infinity],Norm[BF,Infinity]];
        numit += 1;
	];

    If[numit==maxNumIt && BuTools`Verbose, Print["Maximum Number of Iterations reached"]];

    (* Shift Technique *)
    If[shift,
        Lm = Lold;
        If[drift < 0,
            Fm = Fold; (* restore original A2 *)
        ,  (* pos recurrent *)
            G = G + ConstantArray[1,{m,1}].uT;
            Bm = Bold;  (* restore original A0*)
		];
	];
    If[BuTools`Verbose,
        resNorm = Norm[G-Bm-(Lm+Fm.G).G, Infinity];
        Print["Final Residual Error for G: ", resNorm]];

    ret = {};
    Do[
        If[M=="G", AppendTo[ret,G]];
        If[M=="R",
            R = Fm.Inverse[II-(Lm+Fm.G)];
            If[BuTools`Verbose,
                resNorm = Norm[R-Fm-R.(Lm+R.Bm), Infinity];
                Print["Final Residual Error for R: ", resNorm]];
            AppendTo[ret,R]];
        If[M=="U",
            U = Lm+Fm.G;
            If[BuTools`Verbose,
                resNorm = Norm[U-Lm-Fm.Inverse[II-U].Bm,Infinity];
                Print["Final Residual Error for U: ", resNorm]];
            If[continuous, U = lamb*(U-II)];
            AppendTo[ret,U]];
	,{M,Characters[matrices]}];
    If[Length[ret]==1, Return[ret[[1]]], Return[ret]];
];


QBDSolve[B_, L_, F_, L0_, prec_:N[10^-14]]:=
Module[{m, II,lamb,Bm,L0m,pi0,nr,R},
   
    m = Dimensions[L0][[1]];
    II = IdentityMatrix[m];
    R = QBDFundamentalMatrices[B, L, F, "R", prec];
    
	(* Convert to discrete time problem, if needed *)
    If[AnyTrue[Diagonal[L],Negative],
        lamb = Max[-Diagonal[L0]];
        Bm = B / lamb;
        L0m = L0 / lamb + II;
    ,
		Bm = B; L0m=L0];
    
	pi0 = DTMCSolve[L0m+R.Bm];
    nr = Total[pi0.Inverse[II-R]];
    pi0 = pi0 / nr;
	Return[{pi0,R}];
];


QBDStationaryDistr[pi0_, R_, K_]:=
Module[{m, pix,pi},
  
    m = Dimensions[R][[1]];
	pi = Table[0,{(K+1)m}];
	pi[[1;;m]] = pi0;
	pix = pi0;
	Do[
		pix = pix.R;
		pi[[k*m+1;;(k+1)*m]] = pix;
	,{k,K}];
	Return[pi];
];


FluidFundamentalMatrices[Fpp_,Fpm_,Fmp_,Fmm_,matrices_:"P",precision_:N[10^-14],maxNumIt_:50,method_:"ADDA"]:=
Module[{F, sp,sm,II,Pf,A1,A1hat,R1,NN,Xmat,A,B,C,D,S1,U,Ghat,numit,Zm,Psi,
gamma1,gamma2,sA,sD,IA,ID,Dginv,Vginv,Wginv,Eg,Fg,Gg,Hg,diff,neg,nfg,eta,ret,resNorm},

	If[method=="CR",
		F = ArrayFlatten[{{Fpp, Fpm},{Fmp, Fmm}}];
		sp = Dimensions[Fpp][[1]];
		sm = Dimensions[Fmm][[1]];
		II =IdentityMatrix[sm+sp];
		Pf = II + F/Max[-F];
		Zm = Table[0,{sm},{sm}];
		A1 = ArrayFlatten[{{Pf[[1;;sp,1;;sp]]/2, Table[0,{sp},{sm}]},{Pf[[sp+1;;,1;;sp]],Zm}}];
		A1hat = A1;
		R1 = IdentityMatrix[sp]/2;
		NN = ArrayFlatten[{{Pf[[1;;sp,sp+1;;]]/2},{Pf[[sp+1;;,sp+1;;]]}}];
		numit = 0;
		While[Min[Norm[R1,Infinity], Norm[NN,Infinity]]>precision && numit<maxNumIt,
			Xmat = Inverse[II-A1];
			B = Xmat.NN;
			S1 = Xmat[[sp+1;;,1;;sp]].R1;
			U = R1.B[[1;;sp,;;]];
			A1 += ArrayFlatten[{{NN.S1,ArrayFlatten[{{U},{Zm}}]}}];
			A1hat[[1;;sp,sp+1;;]] +=U;
			NN = NN.B[[sp+1;;,;;]];
			R1 = R1.Xmat[[1;;sp,1;;sp]].R1;
			numit++;
		];
		Ghat = Inverse[II-A1hat].ArrayFlatten[{{Pf[[1;;sp,sp+1;;]]/2},{Pf[[sp+1;;,sp+1;;]]}}];
		Psi = Ghat[[1;;sp,;;]];
	, If [method=="ADDA" || method=="SDA",
		A = -Fpp;
        B = Fpm;
        C = Fmp;
        D = -Fmm;
        gamma1 = Max[Diagonal[A]];
        gamma2 = Max[Diagonal[D]];
        If[method=="SDA",
            gamma1 = Max[gamma1,gamma2];
            gamma2 = gamma1;
        ];
        sA = Dimensions[A][[1]];
        sD = Dimensions[D][[1]];
        IA = IdentityMatrix[sA];
        ID = IdentityMatrix[sD];
        A = A + gamma2 IA;
        D = D + gamma1 ID;
        Dginv = Inverse[D];
        Vginv = Inverse[D-C.Inverse[A].B];
        Wginv = Inverse[A-B.Dginv.C];
        Eg = ID - (gamma1+gamma2) Vginv;
        Fg = IA - (gamma1+gamma2) Wginv;
        Gg = (gamma1+gamma2) Dginv.C.Wginv;
        Hg = (gamma1+gamma2) Wginv.B.Dginv;
        diff = 1;
		numit = 0;
        While[diff>precision && numit<maxNumIt,
            Vginv = Eg.Inverse[ID-Gg.Hg];
            Wginv = Fg.Inverse[IA-Hg.Gg];
            Gg += Vginv.Gg.Fg;
            Hg += Wginv.Hg.Eg;
            Eg = Vginv.Eg;
            Fg = Wginv.Fg;
            neg = Norm[Eg,1];
            nfg = Norm[Fg,1];
            If[method=="ADDA",
                eta = Sqrt[nfg/neg];
                Eg = Eg eta;
                Fg = Fg/eta;
                diff = neg nfg;
            ,
                diff = Min[neg,nfg];
            ];
            numit++;
        ];
        Psi = Hg;
	];
	];

    If[numit==maxNumIt && BuTools`Verbose, Print["Maximum Number of Iterations reached"]];

    If[BuTools`Verbose,
        resNorm = Norm[Fpm+Fpp.Psi+Psi.Fmm+Psi.Fmp.Psi, Infinity];
        Print["Final Residual Error for Psi: ", resNorm]];

    ret = {};
    Do[
        If[M=="P", AppendTo[ret,Psi]];
        If[M=="K", AppendTo[ret,Fpp+Psi.Fmp]];
        If[M=="U", AppendTo[ret,Fmm+Fmp.Psi]];
	,{M,Characters[matrices]}];
    If[Length[ret]==1, Return[ret[[1]]], Return[ret]];
];


FluidSolve[Fpp_, Fpm_, Fmp_, Fmm_, prec_:N[10^-14]]:=
Module[{Psi,K,U,mass0,nr,ini,clo},
    {Psi, K, U} = FluidFundamentalMatrices[Fpp, Fpm, Fmp, Fmm, "PKU", prec];
    mass0 = CTMCSolve[U];
    nr = Total[mass0] + 2 Total[mass0.Fmp.Inverse[-K]];
    mass0 = mass0/nr;       
    ini = mass0.Fmp;
    clo = Join[IdentityMatrix[Dimensions[Fpp][[1]]], Psi,2];
	mass0 = Join[Table[0,{Length[Fpp]}],mass0];
	Return[{mass0,ini,K,clo}];
];


FluidStationaryDistr[mass0_, ini_, K_, clo_, x_]:=
Module[{m,y,closing,Ik},
    m = Dimensions[clo][[2]];
    y = Table[0,{Length[x]},{m}];
    closing = Inverse[-K].clo;
	Ik = IdentityMatrix[Dimensions[K][[1]]];
    Do[
        y[[i,All]] = mass0 + ini.(Ik-MatrixExp[K x[[i]]]).closing;
	,{i,Length[x]}];
	Return[y];
];


GeneralFluidSolve[Q_,R_,Q0_:{},prec_:N[10^-14]]:=
Module[{NN,ixz,ixp,ixn,i,Nz,Np,Nn,P,PI,Qv,Rv,
iQv00,Qbar,absRi,Qz,Psi,K,U,Pm,iCn,iCp,clo,
Ua,b,pm,mass0,Qm,Rm,ini,X,Q0v,M,Ma},
	NN=Dimensions[Q][[1]];
	(*Partition the state space according to positive and negative fluid rates*)
	ixz={};
	ixp={};
	ixn={};
	Do[
		If[Abs[R[[i,i]]]<=prec, AppendTo[ixz,i],
			If[R[[i,i]]>prec, AppendTo[ixp,i],
				If[R[[i,i]]<prec, AppendTo[ixn,i]]]];
	,{i,NN}];
	Nz = Length[ixz];
	Np = Length[ixp];
	Nn = Length[ixn];
	(*Create permutation matrix that converts between the original and the partitioned state ordering*)
	P = Table[0,{NN},{NN}];
	Do[P[[i,ixz[[i]]]]=1,{i,Nz}];
	Do[P[[Nz+i,ixp[[i]]]]=1,{i,Np}];
	Do[P[[Nz+Np+i,ixn[[i]]]]=1,{i,Nn}];
	PI=Inverse[P];
	(* reorder states with permutation matrix P: 0, +, - *)
	Qv=P.Q.PI;
	Rv=P.R.PI;
	(* new fluid process censored to states + and - *)
    iQv00 = PseudoInverse[-Qv[[1;;Nz,1;;Nz]]];
    Qbar = Qv[[Nz+1;;, Nz+1;;]] + Qv[[Nz+1;;,1;;Nz]].iQv00.Qv[[1;;Nz,Nz+1;;]];   
    absRi = DiagonalMatrix[Abs[1/Diagonal[Rv[[Nz+1;;,Nz+1;;]]]]];
    Qz = absRi.Qbar;
	(*Calculate fundamental matrices*)
	{Psi,K,U}=FluidFundamentalMatrices[Qz[[1;;Np,1;;Np]],Qz[[1;;Np,Np+1;;]],Qz[[Np+1;;,1;;Np]],Qz[[Np+1;;,Np+1;;]],"PKU",prec];
	(* closing matrix *)
    Pm = Join[IdentityMatrix[Np],Psi,2];
    iCn = absRi[[Np+1;;,Np+1;;]];
    iCp = absRi[[1;;Np,1;;Np]];
    clo = Join[(iCp.Qv[[Nz+1;;Nz+Np,1;;Nz]]+Psi.iCn.Qv[[Nz+Np+1;;,1;;Nz]]).iQv00, Pm.absRi, 2];
	If[Length[Q0]==0, (* regular boundary behavior *)
		clo = clo.P;  (* go back the the original state ordering *)
        (* calculate boundary vector *)
        Ua = iCn.Qv[[Nz+Np+1;;,1;;Nz]].iQv00.Table[1,{Nz},{1}]+ iCn.Table[1,{Nn},{1}] + Qz[[Np+1;;,1;;Np]].Inverse[-K].clo.Table[1,{NN},{1}];
        X = Join[U,Ua,2];
		b = ConstantArray[0,Nn];
		b = Append[b,1];
		{Qm,Rm}=QRDecomposition[Transpose[X]];
		pm=Flatten[Inverse[Rm].Qm.b];
		(* create the result *)
        mass0 = Join[pm.iCn.Qv[[Nz+Np+1;;,1;;Nz]].iQv00, Table[0,{Np}], pm.iCn].P;
        ini = pm.Qz[[Np+1;;,1;;Np]];
	,
	    (* solve a linear system for ini(+), pm(-) and pm(0) *)
        Q0v = P.Q0.PI;
        M = Join[-clo.Rv, Q0v[[Nz+Np+1;;,All]], Q0v[[1;;Nz,All]]];
        Ma = Join[ Inverse[-K].clo.Table[1,{NN},{1}], Table[1,{Nz+Nn},{1}]];
		X = Join[M,Ma,2];
		b = ConstantArray[0,NN];
		b = Append[b,1];
		{Qm,Rm}=QRDecomposition[Transpose[X]];
		pm=Flatten[Inverse[Rm].Qm.b];
        ini = pm[[1;;Np]];
        clo = clo.P;
        mass0 = Join[pm[[Np+Nn+1;;]], Table[0,{Np}], pm[[Np+1;;Np+Nn]]].P;
	];
	Return[{mass0, ini, K, clo}];
];


GM1TypeCaudal[A_,prec_:N[10^-14]]:=
Module[{m,dega,etamin,etamax,eta,i,temp,Dx,Vx,neweta},

	m=Dimensions[A][[1]];
	dega=Dimensions[A][[2]]/m-1;
	etamin=0;
	etamax=1;
	eta=1/2;
	While[etamax - etamin > prec,
		temp = A[[All,dega*m+1;;]];
		Do[temp = temp eta + A[[All,i*m+1;;(i+1)*m]],{i,dega-1,0,-1}];
		neweta=Max[Eigenvalues[temp]];
		If[neweta>eta, etamin=eta, etamax=eta];
        eta=(etamin+etamax)/2;
	];
    {Dx,Vx}=Eigensystem[temp];
	Return[{Re[eta],Re[Vx[[Ordering[Dx][[-1]]]]]}];
];


MG1TypeDecay[A_,prec_:N[10^-14]]:=
Module[{m,dega,etamin,etamax,eta,i,temp,Dx,Vx,neweta},
	
	m=Dimensions[A][[1]];
	dega=Dimensions[A][[2]]/m-1;
	eta=1;
	neweta=0;
	While[neweta - eta < 0,
		eta++;
		temp=A[[All,dega*m+1;;]];
		Do[temp = temp eta + A[[All,i*m+1;;(i+1)*m]],{i,dega-1,0,-1}];
	    neweta=Max[Eigenvalues[temp]];
	];

	etamin=eta-1;
	etamax=eta;
	eta=etamin+1/2;
	While[etamax - etamin > prec,
		temp=A[[All,dega*m+1;;]];
		Do[temp = temp eta + A[[All,i*m+1;;(i+1)*m]],{i,dega-1,0,-1}];
	    neweta=Max[Eigenvalues[temp]];
		If[neweta<eta, etamin=eta, etamax=eta];
        eta=(etamin+etamax)/2;
	];
    {Dx,Vx}=Eigensystem[Transpose[temp]];
	Return[{Re[eta],Re[Vx[[Ordering[Dx][[-1]]]]]}];
];


MG1TypeShifts[A_, shiftType_, prec_:N[10^-14]]:=
Module[{Am,m,v,tau,maxd,II,sumA,beta,i,theta,drift,hatA,uT,rowhatA,colhatA},
	
	m = Dimensions[A][[1]];
	II = IdentityMatrix[m];
	v = Table[0,{m}]; (* default value if redundant *)
	tau = 1; (* default value if redundant *)
	maxd = Dimensions[A][[2]]/m-1;
	sumA = A[[All,maxd*m+1;;]];
	beta = Total[sumA,{2}];
	(* beta = (A_maxd)e + (A_maxd + A_maxd-1)e + ... + (Amaxd+...+A1)e *)
	Do[sumA+=A[[All,i*m+1;;(i+1)*m]];beta+=Total[sumA,{2}],{i,maxd-1,1,-1}];
	sumA += A[[All,1;;m]];
	theta = DTMCSolve[sumA];
	drift = theta.beta;
	(* if drift < 1 : positive recurrent *)
	hatA = Table[0,{m},{m*(maxd+1)}];
	Am = A;
	If[drift < 1,
		If[shiftType=="tau" || shiftType=="dbl",  (* shift tau to infinity *)
			{tau,uT} = MG1TypeDecay[Am,prec];
            Am[[All,m+1;;2*m]] -= II;
            uT /= Total[uT];
            rowhatA = Table[0,{m*(maxd+1)}];
            rowhatA[[maxd*m+1;;]] = uT.Am[[All,maxd*m+1;;]];
			Do[
				rowhatA[[i*m+1;;(i+1)*m]]=tau*rowhatA[[(i+1)*m+1;;(i+2)*m]]+uT.Am[[All,i*m+1;;(i+1)*m]]; 
			,{i,maxd-1,0,-1}];
            hatA = Am - Table[1,{m},{1}].{rowhatA};
		];
		If[shiftType=="dbl", (* shift one to zero *)
			Am=hatA
			(* e is also the right eigenvector of hatA(1) *)
			(* as the shift-tau does not influence G and Ge=e, *)
			(* implying that hatA(1)e = e as G=hatA(G) *)
		];
	    If[shiftType=="one", (* shift one to zero *)
            Am[[All,m+1;;2*m]] -= II;
		];
        If[shiftType=="one" || shiftType=="dbl", (* shift one ot zero *)
            colhatA = Table[0,{m},{maxd+1}];
            colhatA[[All,1]] = Total[Am[[All,1;;m]],{2}]; (* colhatA(:,1) = (A0)e *)
            Do[colhatA[[All,i+1]]=colhatA[[All,i]]+Total[Am[[All,i*m+1;;(i+1)*m]],{2}],{i,maxd}];  (* colhatA(:,i+1) = (A0+A1+...+Ai)e *)
			hatA = Am - KroneckerProduct[colhatA, Table[1,{1},{m}]/m]; (* hatAi = Ai - (A0+A1+...+Ai)e*uT *)
		];
	,
		If[shiftType=="one" || shiftType=="dbl", (* shift one to infinity *)
			Am[[All,m+1;;2*m]] -= II;
			rowhatA = Table[0,{m*(maxd+1)}];
			rowhatA[[maxd*m+1;;]] = theta.Am[[All,maxd*m+1;;]];
			Do[
				rowhatA[[i*m+1;;(i+1)*m]]=rowhatA[[(i+1)*m+1;;(i+2)*m]]+theta.Am[[All,i*m+1;;(i+1)*m]]; 
			,{i,maxd-1,0,-1}];
            hatA = Am - Table[1,{m},{1}].{rowhatA};
		];
		If[shiftType=="dbl", (* shift one to infinity *)
			Am=hatA;
			Am[[All,m+1;;2*m]] = Am[[All,m+1;;2*m]] + II;
			(* v is also the right eigenvector of hatA(tau) *)
			(* as the shift-one does not influence G and Gv=tau*v, *)
			(* implying that hatA(tau)v = tau*v as G=hatA(G) *)
		];
		If[shiftType=="tau" || shiftType=="dbl",  (* shift tau to zero *)
			{tau,v} = GM1TypeCaudal[A,prec];
			Am[[All,m+1;;2*m]] -= II;
			v /= Total[v];
			colhatA = Table[0,{m},{maxd+1}];
			colhatA[[All,1]] = Am[[All,1;;m]].v; (* colhatA(:,1) = (A0)v *)
			Do[colhatA[[All,i+1]]=colhatA[[All,i]]/tau + Am[[All,i*m+1;;(i+1)*m]].v,{i,maxd}];
            hatA = Am - KroneckerProduct[colhatA, Table[1,{1},{m}]];
		];
	];
    hatA[[All,m+1;;2*m]] = hatA[[All,m+1;;2*m]] + II;
	Return[{hatA,drift,tau,v}];
];


MG1FundamentalMatrix[A_, prec_:N[10^-14], maxNumIt_:50, method_:"ShiftPWCR", maxNumRoot_:2048, shiftType_:"one"]:=
Module[{mode,Dm,Dold,drift,tau,v,lastiter,m,G,Aeven,Aodd,
Ahatodd,Ahateven,Rj,numit,nj,temp,temp1,temp2,temp3,temp4,
Ahatnew,Anew,Aodd1,nAnew,nAhatnew,deg,deghat,stopv,
vec1,vec2,Rnewj,Gold,resNorm},

	If[method=="ShiftPWCR", mode="ShiftPWCR", mode="PWCR"];
	If[method!="CR" && method!="ShiftPWCR" && BuTools`Verbose,
		Print["MG1FundamentalMatrix: method not supported. Falling back to ShiftPWCR."];
		mode="ShiftPWCR"];

	If[Length[Dimensions[A]]>2,Dm = ArrayFlatten[{A}],Dm=A];

	If[mode=="ShiftPWCR",
		If[BuTools`Verbose, Dold = Dm];
		{Dm,drift,tau,v} = MG1TypeShifts[Dm, shiftType, prec];
	];

	(* start Cyclic Reduction *)
	lastiter = 0;
	m = Dimensions[Dm][[1]];
	Dm = Transpose[Dm];
	Dm = Join[Dm, Table[0, {(2^(1+Floor[Log2[Dimensions[Dm][[1]]/m-1]])+1)*m-Dimensions[Dm][[1]]},{m}]];

	(* Step 0 *)
	G = Table[0,{m},{m}];

	Aeven = Pick[Dm,Mod[KroneckerProduct[{Range[Dimensions[Dm][[1]]/m]},Table[1,{m}]][[1]],2],1];
	Aodd = Pick[Dm,Mod[KroneckerProduct[{Range[Dimensions[Dm][[1]]/m]},Table[1,{m}]][[1]],2],0];
	Ahatodd = Join[Aeven[[m+1;;,All]], Dm[[-m;;,All]]];
	Ahateven = Aodd;

	Rj = Dm[[m+1;;2*m,All]];
	Do[Rj+=Dm[[(i-1)*m+1;;i*m,All]],{i,3,Dimensions[Dm][[1]]/m}];
	Rj = Inverse[IdentityMatrix[m]-Rj];
	Rj = Dm[[1;;m,All]].Rj;

	numit=0;
	While[numit<maxNumIt,
		numit++;
		nj = Dimensions[Aodd][[1]]/m-1;
		If[nj>0,
			(* Evaluate the 4 functions in the nj+1 roots using FFT *)
			(* prepare for FFTs (such that they can be performed in 4 calls) *)
			temp1 = ArrayReshape[Aodd[[1;;(nj+1)*m,All]],{nj+1,m^2}];
			temp2 = ArrayReshape[Aeven[[1;;(nj+1)*m,All]],{nj+1,m^2}];
			temp3 = ArrayReshape[Ahatodd[[1;;(nj+1)*m,All]],{nj+1,m^2}];
			temp4 = ArrayReshape[Ahateven[[1;;(nj+1)*m,All]],{nj+1,m^2}];
			(* FFTs *)
			temp1 = Transpose[Map[Fourier[#, FourierParameters -> {1, -1}]&,Transpose[temp1]]];
			temp2 = Transpose[Map[Fourier[#, FourierParameters -> {1, -1}]&,Transpose[temp2]]];
			temp3 = Transpose[Map[Fourier[#, FourierParameters -> {1, -1}]&,Transpose[temp3]]];
			temp4 = Transpose[Map[Fourier[#, FourierParameters -> {1, -1}]&,Transpose[temp4]]];
			(* reform the 4*(nj+1) matrices *)
			temp1 = ArrayReshape[temp1,{m*(nj+1),m}];
			temp2 = ArrayReshape[temp2,{m*(nj+1),m}];
			temp3 = ArrayReshape[temp3,{m*(nj+1),m}];
			temp4 = ArrayReshape[temp4,{m*(nj+1),m}];
			(* Next, we perform a point-wise evaluation of (6.20) - Thesis Meini *)
			Ahatnew = Table[0,{(nj+1)*m},{m}];
			Anew = Table[0,{(nj+1)*m},{m}];
			Do[
				Ahatnew[[(cnt-1)*m+1;;cnt*m,1;;m]] = temp4[[(cnt-1)*m+1;;cnt*m,All]] + temp2[[(cnt-1)*m+1;;cnt*m,All]].Inverse[IdentityMatrix[m]-temp1[[(cnt-1)*m+1;;cnt*m,All]]].temp3[[(cnt-1)*m+1;;cnt*m,All]];
				Anew[[(cnt-1)*m+1;;cnt*m,1;;m]] = Exp[-(cnt-1)*2*I*Pi/(nj+1)] temp1[[(cnt-1)*m+1;;cnt*m,All]] + temp2[[(cnt-1)*m+1;;cnt*m,All]].Inverse[IdentityMatrix[m]-temp1[[(cnt-1)*m+1;;cnt*m,All]]].temp2[[(cnt-1)*m+1;;cnt*m,All]];
			,{cnt,nj+1}];
			(* We now invert the FFTs to get Pz and Phatz *)
			(* prepare for IFFTs (in 2 calls) *)
			Ahatnew = ArrayReshape[Ahatnew[[1;;(nj+1)*m,All]],{(nj+1),m^2}];
			Anew = ArrayReshape[Anew[[1;;(nj+1)*m,All]],{(nj+1),m^2}];
			(* IFFTs *)
			Ahatnew=Re[Transpose[Map[InverseFourier[#, FourierParameters -> {1, -1}]&,Transpose[Ahatnew]]]];
			Anew=Re[Transpose[Map[InverseFourier[#, FourierParameters -> {1, -1}]&,Transpose[Anew]]]];
			(* reform matrices Pi and Phati *)
			Ahatnew = ArrayReshape[Ahatnew,{m*(nj+1),m}];
			Anew = ArrayReshape[Anew,{m*(nj+1),m}];
		,
			(* series Aeven, Aodd, Ahateven and Ahatodd are constant *)
			temp = Aeven.Inverse[IdentityMatrix[m]-Aodd];
			Ahatnew = Ahateven + temp.Ahatodd;
			Anew = Join[temp.Aeven, Aodd];
			Aodd1 = Aodd;
		];

		nAnew = 0;
		deg = Dimensions[Anew][[1]]/m;
		Do[nAnew = Max[nAnew,Norm[Anew[[i*m+1;;(i+1)*m,All]],Infinity]],{i,deg/2,deg-1}];
		nAhatnew = 0;
		deghat = Dimensions[Ahatnew][[1]]/m;
		Do[nAhatnew = Max[nAhatnew,Norm[Ahatnew[[i*m+1;;(i+1)*m,All]],Infinity]],{i,deghat/2,deghat-1}];
    
		(* c) the test *)
		While[ (nAnew>(nj+1)*prec || nAhatnew>(nj+1)*prec) && nj+1 < maxNumRoot,
			nj = 2*(nj+1)-1;
			stopv = Min[nj+1, Dimensions[Aodd][[1]]/m];

			(* prepare for FFTs *)
			temp1 = ArrayReshape[Aodd[[1;;stopv*m,All]],{stopv,m^2}];
			temp2 = ArrayReshape[Aeven[[1;;stopv*m,All]],{stopv,m^2}];
			temp3 = ArrayReshape[Ahatodd[[1;;stopv*m,All]],{stopv,m^2}];
			temp4 = ArrayReshape[Ahateven[[1;;stopv*m,All]],{stopv,m^2}];
			(* FFTs *)
			temp1 = Transpose[Map[Fourier[PadRight[#,nj+1,0], FourierParameters -> {1, -1}]&,Transpose[temp1]]];
			temp2 = Transpose[Map[Fourier[PadRight[#,nj+1,0], FourierParameters -> {1, -1}]&,Transpose[temp2]]];
			temp3 = Transpose[Map[Fourier[PadRight[#,nj+1,0], FourierParameters -> {1, -1}]&,Transpose[temp3]]];
			temp4 = Transpose[Map[Fourier[PadRight[#,nj+1,0], FourierParameters -> {1, -1}]&,Transpose[temp4]]];
			(* reform the 4*(nj+1) matrices *)
			temp1 = ArrayReshape[temp1,{m*(nj+1),m}];
			temp2 = ArrayReshape[temp2,{m*(nj+1),m}];
			temp3 = ArrayReshape[temp3,{m*(nj+1),m}];
			temp4 = ArrayReshape[temp4,{m*(nj+1),m}];
			(* Next, we perform a point-wise evaluation of (6.20) - Thesis Meini *)
			Ahatnew = Table[0,{(nj+1)*m},{m}];
			Anew = Table[0,{(nj+1)*m},{m}];
			Do[
				Ahatnew[[(cnt-1)*m+1;;cnt*m,1;;m]] = temp4[[(cnt-1)*m+1;;cnt*m,All]] + temp2[[(cnt-1)*m+1;;cnt*m,All]].Inverse[IdentityMatrix[m]-temp1[[(cnt-1)*m+1;;cnt*m,All]]].temp3[[(cnt-1)*m+1;;cnt*m,All]];
				Anew[[(cnt-1)*m+1;;cnt*m,1;;m]] = Exp[-(cnt-1)*2*I*Pi/(nj+1)] temp1[[(cnt-1)*m+1;;cnt*m,All]] + temp2[[(cnt-1)*m+1;;cnt*m,All]].Inverse[IdentityMatrix[m]-temp1[[(cnt-1)*m+1;;cnt*m,All]]].temp2[[(cnt-1)*m+1;;cnt*m,All]];
			,{cnt,nj+1}];
			(* We now invert the FFTs to get Pz and Phatz *)
			(* prepare for IFFTs (in 2 calls) *)
			Ahatnew = ArrayReshape[Ahatnew[[1;;(nj+1)*m,All]],{(nj+1),m^2}];
			Anew = ArrayReshape[Anew[[1;;(nj+1)*m,All]],{(nj+1),m^2}];
			(* IFFTs *)
			Ahatnew=Re[Transpose[Map[InverseFourier[#, FourierParameters -> {1, -1}]&,Transpose[Ahatnew]]]];
			Anew=Re[Transpose[Map[InverseFourier[#, FourierParameters -> {1, -1}]&,Transpose[Anew]]]];
			(* reform matrices Pi and Phati *)
			Ahatnew = ArrayReshape[Ahatnew,{m*(nj+1),m}];
			Anew = ArrayReshape[Anew,{m*(nj+1),m}];
	
			vec1 = Table[0,{m}];
			vec2 = Table[0,{m}];
			Do[
				vec1 += i*Total[Anew[[i*m+1;;(i+1)*m,All]],{1}];
				vec2 += i*Total[Ahatnew[[i*m+1;;(i+1)*m,All]],{1}];
			,{i,Dimensions[Anew][[1]]/m-1}];
			nAnew = 0;
			deg = Dimensions[Anew][[1]]/m;
			Do[nAnew = Max[nAnew,Norm[Anew[[i*m+1;;(i+1)*m,All]],Infinity]],{i,deg/2,deg-1}];
			nAhatnew = 0;
			deghat = Dimensions[Ahatnew][[1]]/m;
			Do[nAhatnew = Max[nAhatnew,Norm[Ahatnew[[i*m+1;;(i+1)*m,All]],Infinity]],{i,deghat/2,deghat-1}];
		];
		If[(nAnew>(nj+1)*prec || nAhatnew>(nj+1)*prec) && nj+1>=maxNumRoot,
			Print["MG1FundamentalMatrix: Maximum number of iterations reached, accuracy might be affected!"]];
		If[nj>1,
			Anew = Anew[[1;;m*(nj+1)/2,All]];
			Ahatnew = Ahatnew[[1;;m*(nj+1)/2,All]];
		];
		(* compute Aodd, Aeven, ... *)
		Aeven = Pick[Anew,Mod[KroneckerProduct[{Range[Dimensions[Anew][[1]]/m]},Table[1,{m}]][[1]],2],1];
		Aodd = Pick[Anew,Mod[KroneckerProduct[{Range[Dimensions[Anew][[1]]/m]},Table[1,{m}]][[1]],2],0];
		Ahateven = Pick[Ahatnew,Mod[KroneckerProduct[{Range[Dimensions[Ahatnew][[1]]/m]},Table[1,{m}]][[1]],2],1];
		Ahatodd = Pick[Ahatnew,Mod[KroneckerProduct[{Range[Dimensions[Ahatnew][[1]]/m]},Table[1,{m}]][[1]],2],0];
		
		If[BuTools`Verbose,Print["The evaluation of the iteration required ",nj+1," roots\n"]];
        
		(* test stopcriteria *)
		If[mode=="PWCR" || mode=="DCR",
			Rnewj = Anew[[m+1;;2*m,All]];
			Do[Rnewj += Anew[[(i-1)*m+1;;i*m,All]],{i,3,Dimensions[Anew][[1]]/m}];
			Rnewj = Inverse[IdentityMatrix[m]-Rnewj];
			Rnewj = Anew[[1;;m,All]].Rnewj;
			If[Max[Abs[Rj-Rnewj]] < prec || Max[Total[IdentityMatrix[m]-Anew[[1;;m,All]].Inverse[IdentityMatrix[m]-Anew[[m+1;;2*m,All]]],{1}]] < prec,
				G = Ahatnew[[1;;m,All]];
				Do[G+=Rnewj.Ahatnew[[(i-1)*m+1;;i*m,All]],{i,2,Dimensions[Ahatnew][[1]]/m}];
				G = Dm[[1;;m,All]].Inverse[IdentityMatrix[m]-G];
				Break[];
			];
			Rj=Rnewj;
			(* second condition tests whether Ahatnew is degree 0 (numerically) *)
			If[Norm[Anew[[1;;m,1;;m]]]<prec || Total[Ahatnew[[m+1;;,All]]]<prec || Max[Total[IdentityMatrix[m]-Dm[[1;;m,All]].Inverse[IdentityMatrix[m]-Ahatnew[[1;;m,All]]],{1}]] < prec,
				G = Dm[[1;;m,All]].Inverse[IdentityMatrix[m]-Ahatnew[[1;;m,All]]];
				Break[];
			];
        ,
			Gold = G;
			G = Dm[[1;;m,All]].Inverse[IdentityMatrix[m]-Ahatnew[[1;;m,All]]];
			If[Norm[G-Gold,Infinity]<prec || Norm[Ahatnew[[m+1;;,All]],Infinity]<prec, Break[]];
		];   
	];
	If[numit==maxNumIt && G BuTools`Verbose, 
		Print["Maximum Number of Iterations reached"];
		G = Dm[[1;;m,All]].Inverse[IdentityMatrix[m]-Ahatnew[[1;;m,All]]];
	];
	G = Transpose[G];
	If[mode=="ShiftPWCR",
		If[shiftType=="one", G+=Boole[drift<1]*Table[1,{m},{m}]/m,
			If[shiftType=="tau", G+=Boole[drift>1]*tau*Transpose[{v}].Table[1,{1},{m}],
				If[shiftType=="dbl", G+=Boole[drift<1]*Table[1,{m},{m}]/m + Boole[drift>1]*tau*Transpose[{v}].Table[1,{1},{m}]]]]];
            
	If[BuTools`Verbose,
		If[mode=="PWCR", Dm=Transpose[Dm], Dm = Dold];
		temp = Dm[[All,-m;;]];
		Do[temp=Dm[[All,(i-1)*m+1;;i*m]]+temp.G,{i,Dimensions[Dm][[2]]/m-1,1,-1}];
		resNorm = Norm[G-temp,Infinity];
		Print["Final Residual Error for G: ",resNorm];
	];
	Return[G];
];


MG1StationaryDistr[A_, B_:{}, G_:{}, K_:500, prec_:N[10^-14]]:=
Module[{Am,m,II,dega,Bm,mb,degb,sumA,beta,theta,drift,Gm,g,
sumBB0,pi0,Bbeta,Km,kappa,psi1,psi2,temp,tildekappa1,
invbarA1,sumpi,numit,pi,pix},

	(* Compute g *)
    If[Length[G]==0, Gm = MG1FundamentalMatrix[A, prec],Gm=G];
    g = DTMCSolve[Gm];
    
	Am = ArrayFlatten[{A}];
    m = Dimensions[Am][[1]];
    II = IdentityMatrix[m];
    dega = Dimensions[Am][[2]]/m-1;
    If[Length[B]==0,
		mb = m;
		degb = dega;
		Bm = Am;
	,
		Bm = ArrayFlatten[{B}];
		mb = Dimensions[Bm][[1]];
		If[mb!=m,Throw["MG1StationaryDistr: Matrix B has an incorrect number of columns"]];
		degb = Floor[(Dimensions[Bm][[2]]-mb)/m];
	];
    
    (* Compute theta and beta, sum_v>=0 Av and sum_v>=k Av G^v-1 *)
    (* the last sums (for k=1,...,amax) are stored in A *)
	sumA = Am[[All,dega*m+1;;]];
	beta = Total[sumA,{2}];
	(* beta = (A_maxd)e + (A_maxd + A_maxd-1)e + ... + (Amaxd+...+A1)e *)
	Do[
		sumA += Am[[All,i*m+1;;(i+1)*m]];
		Am[[All,i*m+1;;(i+1)*m]] += Am[[All,(i+1)*m+1;;(i+2)*m]].G;
		beta += Total[sumA,{2}];
	,{i,dega-1,1,-1}];
	sumA += Am[[All,1;;m]];
	theta = DTMCSolve[sumA];
	drift = theta.beta;
	If[drift >= 1,Throw["MG1StationaryDistr: the Markov chain characterized by A is not positive recurrent"]];
    If[Length[Bm]==0, 
		pi0=(1-drift)*g;
	,
		(* Compute sum_v>=1 Bv, sum_v>=1 (v-1) Bv e, sum_v>=k Bv G^v-1 *)
		(* the last sums (for k=1,...,bmax) are stored in B *)
        sumBB0 = Bm[[All,mb+(degb-1)*m+1;;]];
        Bbeta = Table[0,{mb}];
		Do[
			Bbeta += Total[sumBB0,{2}];
            sumBB0 += Bm[[All,mb+(i-1)*m+1;;mb+i*m]];
            Bm[[All,mb+(i-1)*m+1;;mb+i*m]] += Bm[[All,mb+i*m+1;;mb+(i+1)*m]].G;
		,{i,degb-1,1,-1}];        
		(* Compute K, kappa *)
        Km = Bm[[All,1;;mb]]+Bm[[All,mb+1;;mb+m]].G;
        kappa = DTMCSolve[Km];
		(* Compute psi1, psi2 *)
        temp = Total[Inverse[II-sumA-Transpose[{Table[1,{m}]-beta}].{g}],{2}];
        psi1=(II-Am[[All,1;;m]]-Am[[All,m+1;;2*m]]).temp + (1-drift)^(-1)*Total[Am[[All,1;;m]],{2}];
        psi2=Table[1,{mb}] + (sumBB0-Bm[[All,mb+1;;mb+m]]).temp+(1-drift)^(-1)*Bbeta;
        (* Compute kappa1 *)
		tildekappa1=psi2 + Bm[[All,mb+1;;mb+m]].Inverse[II-Am[[All,m+1;;2*m]]].psi1;
        (* Compute pi_0 *)
        pi0=(kappa.tildekappa1)^(-1)*kappa;
	];
	(* Start stable RAMASWAMI formula *)
	numit = 1;
	sumpi = Total[pi0];
	pi = {};
	invbarA1 = Inverse[II-Am[[All,m+1;;2*m]]];
	While[sumpi< 1-10^-10 && numit < K,
		If[numit <= degb,
			If[Length[Bm]==0, 
				pix = pi0.Am[[All,numit*mb+1;;(numit+1)*mb]];
			,
				pix = pi0.Bm[[All,mb+(numit-1)*m+1;;mb+numit*m]];
			];
		,
			pix = Table[0,{m}];
		];
		Do[pix+=pi[[numit-j]].Am[[All,(j+1)*m+1;;(j+2)*m]],{j,1,Min[numit-1,dega-1]}];
		pix = pix.invbarA1;
		sumpi += Total[pix];
		AppendTo[pi,pix];
		If[BuTools`Verbose, Print["Accumulated mass of the first ",numit, " (reblocked) components:", sumpi]];
        numit++;
	];
	PrependTo[pi,pi0];
	If[numit == K && BuTools`Verbose, Print["Maximum Number of Components reached"]];
    Return[Flatten[pi]];
];


GM1FundamentalMatrix[A_, prec_:N[10^-14], maxNumIt_:50, method_:"ShiftPWCR", dual_:"R", maxNumRoot_:2048, shiftType_:"one"]:=
Module[{m,dega,II,sumA,beta,theta,drift,Am,R,G,eta,v,sumAeta,ci},

	Am = ArrayFlatten[{A}];
	m = Dimensions[Am][[1]];
	II = IdentityMatrix[m];
	dega = Dimensions[Am][[2]]/m-1;

	(* compute invariant vector of A and the drift *)
	(* drift > 1: positive recurrent GIM1, drift < 1: transient GIM1 *)
	sumA = Am[[All,dega*m+1;;]];
	beta = Total[sumA,{2}];
	(* beta = (A_maxd)e + (A_maxd + A_maxd-1)e + ... + (Amaxd+...+A1)e *)
	Do[
		sumA += Am[[All,i*m+1;;(i+1)*m]];
		beta += Total[sumA,{2}];
	,{i,dega-1,1,-1}];
	sumA += Am[[All,1;;m]];
	theta = DTMCSolve[sumA];
	drift = theta.beta;
	If[dual=="R" || (dual=="A" && drift<=1), (* RAM dual *)
		(* compute the RAM Dual process *)
		Do[Am[[All,i*m+1;;(i+1)*m]]=DiagonalMatrix[1/theta].Transpose[Am[[All,i*m+1;;(i+1)*m]]].DiagonalMatrix[theta],{i,0,dega}];
	,
		(* Bright dual *)
		If[drift > 1, (* A -> positive recurrent GIM1 *)
			(* compute the Caudal characteristic of A *)
			{eta,v}=GM1TypeCaudal[Am];
		,
			(* A -> transient GIM1 (=recurrent MG1) *)
			{eta,v}=MG1TypeDecay[Am];
		];
		(* compute invariant vector of A0+A1*eta+A2*eta^2+...+Amax*eta^max *)
		sumAeta = eta^dega * Am[[All,dega*m+1;;]];
		Do[sumAeta+=eta^i * Am[[All,i*m+1;;(i+1)*m]],{i,dega-1,0,-1}];
		ci = BuTools`CheckInput;
		BuTools`CheckInput = False;
		theta = DRPSolve[sumAeta+(1-eta)*II];
		BuTools`CheckInput = ci;
		(* compute the Bright Dual process *)
		Do[Am[[All,i*m+1;;(i+1)*m]]=eta^(i-1) * DiagonalMatrix[1/theta].Transpose[Am[[All,i*m+1;;(i+1)*m]]].DiagonalMatrix[theta],{i,0,dega}];
	];
    G = MG1FundamentalMatrix[Am, prec, maxNumIt, method, maxNumRoot, shiftType];

	If[dual=="R" || (dual=="A" && drift<=1), (* RAM dual *)
		R = DiagonalMatrix[1/theta].Transpose[G].DiagonalMatrix[theta];
	, (* Bright dual *)
		R = DiagonalMatrix[1/theta].Transpose[G].DiagonalMatrix[theta]*eta;
	];   
	Return[R];
];


GM1StationaryDistr[B_, R_, K_]:=
Module[{m,II,temp,maxb,BR,pi,sumpi,numit,pix},

    m = Dimensions[R][[1]];
    II = IdentityMatrix[m];
    temp = Inverse[II-R];
	If[AnyTrue[temp,#<-100*BuTools`CheckPrecision &], Throw["GM1StationaryDistr: The spectral radius of R is not below 1: the Markov chain is not pos. recurrent"]];
	maxb = Length[B];
    BR=B[[maxb]];
	Do[BR=R.BR+B[[i]],{i,maxb-1,1,-1}];
    pix = DTMCSolve[BR]; (* compute pi_0 *)
    pix /= Total[pix.temp]; (* normalize pi_0 *)
    sumpi = Total[pix];
    numit = 1;
	pi = {pix};
	While[sumpi<1-10^-10 && numit<1+K,
		pix = pix.R; (* compute pi_(numit+1) *)
        numit++;
        sumpi += Total[pix];
		AppendTo[pi,pix];
        If[BuTools`Verbose, Print["Accumulated mass after ",numit," iterations: ",sumpi]];
	];
    If[BuTools`Verbose && numit == K+1, Print["Maximum Number of Components ", numit-1, " reached"]];
	Return[Flatten[pi]];
];


End[(* Private *)];
EndPackage[];
