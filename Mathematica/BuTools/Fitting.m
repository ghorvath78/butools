(* ::Package:: *)

(*
   BuTools Fitting Package
*)

BeginPackage["BuTools`Fitting`"];
SquaredDifference::usage = "sd = SquaredDifference[p1, p2]: Returns the squared difference between two vectors";
RelativeEntropy::usage = "re = RelativeEntropy[p1, p2]: Returns the relative entropy (aka Kullback–Leibler divergence) of two vectors";
EmpiricalSquaredDifference::usage = "sd = EmpiricalSquaredDifference[f1, f2, intBounds]: Returns the squared difference of two continuous functions given by samples and the bounds of the corresponding intervalls";
EmpiricalRelativeEntropy::usage = "re = EmpiricalRelativeEntropy[f1, f2, intBounds]: Returns the relative entropy (aka Kullback–Leibler divergence) of two continuous functions given by samples and the bounds of the corresponding intervalls";
LikelihoodFromTrace::usage = "logli = LikelihoodFromTrace[trace, X, Y, prec]: Evaluates the log-likelihood of a trace with the given PH distribution or MAP";
PHFromTrace::usage = "{alpha, A, logli} = PHFromTrace[trace, orders, maxIter, stopCond, initial, result]: Performs PH distribution fitting using the EM algorithm";
MAPFromTrace::usage = "{D0, D1, logli} = MAPFromTrace[trace, orders, maxIter, stopCond, initial, result]: Performs MAP fitting using the EM algorithm";


Begin["`Private`"];


Needs["BuTools`MC`"];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckInput"]],BuTools`CheckInput=True];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`CheckPrecision"]],BuTools`CheckPrecision=N[10^-12]];
If[Not[MemberQ[Names["BuTools`*"],"BuTools`Verbose"]],BuTools`Verbose=False];


SquaredDifference[p1_, p2_]:=Dot[p1-p2,p1-p2];


RelativeEntropy[p1_, p2_]:=
Module[{re},
	re = 0;
	Do[If[p1[[i]]>0,re+=p1[[i]]*Abs[Log[p1[[i]]/p2[[i]]]]],{i,Length[p1]}];
    Return[re];
];


EmpiricalSquaredDifference[f1_, f2_, intBounds_]:=
Module[{intlens,p1,p2},
    intlens = intBounds[[2;;]] - intBounds[[;;-2]];
    p1 = f1 * intlens;
    p2 = f2 * intlens;
    Return[SquaredDifference[p1, p2]];
];


EmpiricalRelativeEntropy[f1_, f2_, intBounds_]:=
Module[{intlens,p1,p2},
    intlens = intBounds[[2;;]] - intBounds[[;;-2]];
    p1 = f1 * intlens;
    p2 = f2 * intlens;
    Return[RelativeEntropy[p1, p2]];
];


LikelihoodFromTrace[trace_,X_,Y_,prec_:N[10^-14]]:=
Module[{alpha, A, trc,lambda,P,a,eps,lpoi,poi,spoi,fx,k,first,coeffv,
maxIter,D0,D1,NN,L,scale,sc,l,logli,ix,ixrev},
	L = Length[trace];
    If[Length[Dimensions[X]]==1,
        (* We have a PH distribution. We can sort it to make the computation faster *)
        alpha = X;
        A = Y;
        trc = Sort[trace];
        lambda = Max[Abs[Diagonal[A]]];
        P = A/lambda + IdentityMatrix[Dimensions[A][[1]]];
        a = Total[-A,{2}];
        eps = Max[prec, 10^(Log10[prec] + Log10[lambda])];
        lpoi = -lambda*trc;
        trc = Log[trc];
        poi =  Exp[lpoi];
        spoi = poi;
        fx = poi*(alpha.a);
        k = 1;
        first = 1;
        coeffv = alpha;
        maxIter = 10000;
        While[first<=L && k<maxIter,
            coeffv = coeffv.P;
            lpoi[[first;;]] += Log[lambda] + trc[[first;;]] - Log[k];
            poi[[first;;]] = Exp[lpoi[[first;;]]];
            spoi[[first;;]] += poi[[first;;]];
            fx[[first;;]] += poi[[first;;]] * coeffv.a;
            k++;
            first += LengthWhile[spoi[[first;;]], #>1-eps&];
        ];
        logli = Total[Log[fx]]/L;
    ,        
        D0 = X;
        D1 = Y;
        NN = Dimensions[D0][[1]];
        (* first we calculate matrix e^(D0*x(i))*D1 for each sample *)
		ix = Ordering[trace];
		trc = trace[[ix]];
        lambda = Max[Abs[Diagonal[D0]]];
        P = D0/lambda + IdentityMatrix[NN];
        eps = Max[prec, 10^(Log10[prec] + Log10[lambda])];
        lpoi = -lambda*trc;
        trc = Log[trc];
        poi =  Exp[lpoi];
        spoi = poi;
        coeffv = D1;
        fx = KroneckerProduct[{poi},coeffv];
        k = 1;
        first = 1;
        maxIter = 10000;
        While[first<=L && k<maxIter,
            coeffv = P.coeffv;
            lpoi[[first;;]] += Log[lambda] + trc[[first;;]] - Log[k];
            poi[[first;;]] = Exp[lpoi[[first;;]]];
            spoi[[first;;]] += poi[[first;;]];
            fx[[All,(first-1)*NN+1;;]] += KroneckerProduct[{poi[[first;;]]},coeffv];
            k++;
			first += LengthWhile[spoi[[first;;]], #>1-eps&];
        ];
        alpha = DTMCSolve[Inverse[-D0].D1];
        l = alpha;
        sc = 0;
        ixrev=Ordering[ix];
        Do[
            l = l.fx[[All,(ixrev[[i]]-1)*NN+1;;ixrev[[i]]*NN]];
            If[Mod[i,10]==0,
                (* sometimes we need to rescale the results to avoid "nan"s *)
                scale = Ceiling[Log2[Total[l]]];            
                If[scale>1, l = l/2^scale; sc += scale];
                If[scale<-10, scale += 10; l = l/2^scale; sc += scale];
            ];
        ,{i,L}];
        logli = (Log[Total[l]]+sc*Log[2]) / L;
    ];
	Return[logli];
];


AllOrders[branches_,sumorders_]:=
Module[{o,x,xt,found},
    If[branches==1,
        o={{sumorders}};
    ,
        o = {};
        Do[x = AllOrders[branches-1,sumorders-ii];
            Do[xt = Sort[Append[x[[j]],ii]];
                (* check if we have it already *)
                found = False;
                Do[If[o[[ok]]==xt, found=True; Break[]],{ok,Length[o]}];
                If[Not[found], AppendTo[o,xt]];
            ,{j,Length[x]}];
        ,{ii,sumorders-branches+1}];
    ];
	Return[o];
];


PHFromTrace[trace_,orders_,maxIter_:200,stopCond_:N[10^-7],initial_:{},result_:"vecmat"]:=
Module[{bestAlpha,bestA,bestLogli,bestOrders,allord,alpha,A,l,
M,K,alphav,lambda,trm,inim,logli,ologli,Q,steps,nor,v1,v2,t1,ix,
NN,blk},
    If [NumberQ[orders],
        bestAlpha = {};
        bestA = {};
        bestLogli = -Infinity;
        bestOrders= {};
		Do[
            allord = AllOrders[br, orders];
            Do[
                {alpha,A,l} = PHFromTrace[trace, allord[[ox]], maxIter, stopCond, initial, result];
                If[l>bestLogli,
                    bestAlpha = alpha;
                    bestA = A;
                    bestLogli = l;
                    bestOrders = allord[[ox]];
                ];
            ,{ox,Length[allord]}];
        ,{br,2,orders}];
        If[BuTools`Verbose,Print["Best solution: logli=",bestLogli,", orders=",bestOrders]];
		alpha = bestAlpha;
		A = bestA;
		logli = bestLogli;
    ,
		M = Length[orders];
		K = Length[trace];
		(* initial alpha and lambda is such that the mean is matched *)
		If[Length[initial]==0,
			alphav = Table[1./M, {M}];
			lambda = DiagonalMatrix[orders].Range[M];
			trm = Mean[trace];
			inim = Total[alphav/Range[M]];
			lambda = lambda * inim / trm;
		,
		If[Length[initial]==2,
			If[Length[initial[[1]]]==M && Length[initial[[2]]]==M,
				alphav = initial[[1]];
				lambda = initial[[2]];
			,
				Throw["The length of the initial branch probability and rate vectors is not consistent with the length of the orders vector!"];
			];
		,
		Throw["Invalid initial branch probability and rate vectors!"];
		]];

		ologli = 1;
		logli = 0.00000000001;
		steps = 1;
		t1 = AbsoluteTime[];
		While[Abs[(ologli-logli)/logli]>stopCond && steps<=maxIter,
			ologli = logli;
			(* E-step: *)
			Q={};
			Do[AppendTo[Q,(alphav[[i]]*(lambda[[i]]*trace)^(orders[[i]]-1) / (orders[[i]]-1)! * lambda[[i]]) * Exp[-lambda[[i]]*trace]],{i,M}];
			logli = Total[Log[Total[Q,{1}]]] / K;
			nor = Total[Q,{1}];
			Do[Q[[i,All]]=Q[[i,All]]/nor,{i,M}];
			(* M-step: *)
			v1 = Total[Q,{2}];
			v2 = Q.trace;
			alphav = v1/K;
			lambda = DiagonalMatrix[orders].v1 / v2;
			steps++;
			If[BuTools`Verbose && AbsoluteTime[]-t1>2,
				Print["Num of iterations: ",steps,", logli: ", logli];
				t1 = AbsoluteTime[];
			];
		];
		If[BuTools`Verbose,
			Print["Num of iterations: ",steps,", logli: ", logli];
			Print["EM algorithm terminated. (orders=",orders,")"];
        ];

		If[result=="vecvec",
			alpha = alphav;
			A = lambda;
		, If[result=="vecmat",
			NN = Total[orders];
			alpha = Table[0,{NN}];
			A = Table[0,{NN},{NN}];
			ix = 1;
			Do[
				alpha[[ix]] = alphav[[i]];
				blk = -DiagonalMatrix[Table[1,{orders[[i]]}]];
				If[orders[[i]]>1, blk+=DiagonalMatrix[Table[1,{orders[[i]]-1}],1]];
				A[[ix;;ix+orders[[i]]-1, ix;;ix+orders[[i]]-1]] = lambda[[i]]*blk;
				ix += orders[[i]];
			,{i,M}];
        ]];
	];
    Return[{alpha,A,logli}];
];


MAPFromTrace[trace_,orders_,maxIter_:200,stopCond_:N[10^-7],initial_:{},result_:"matmat"]:=
Module[{bestX,bestY,bestLogli,bestOrders,allord,oX,oY,
Ascale,Bscale,logli,D0,D1,steps,P,A,B,K,M,alphav,lambda,trm,inim,
ologli,t1,prev,scprev,scale,Av,Ascalev,Bv,Bscalev,next,llh,illh,
AB,v1,v2,Avv,nor,NN,ix,blk,indicesTo,indicesFrom,Q,l},

    If[NumberQ[orders],
        bestX = {};
        bestY = {};
        bestLogli = -Infinity;
        bestOrders= {};
        Do[
            allord = AllOrders[br,orders];
            Do[
                If[BuTools`Verbose, Print["Trying orders ",allord[[ox]],"..."]];
                {oX,oY,l} = MAPFromTrace[trace, allord[[ox]], maxIter, stopCond, initial, result];
                If[l > bestLogli,
                    bestX = oX;
                    bestY = oY;
                    bestLogli = l;
                    bestOrders = allord[[ox]];
                ];
            ,{ox,Length[allord]}];
        ,{br,2,orders}];
        D0 = bestX;
        D1 = bestY;
        logli = bestLogli;
        If[BuTools`Verbose, Print["Best solution: logli=",bestLogli,", orders=",bestOrders]];
	,
		M = Length[orders];
		K = Length[trace];
		(* initial alpha and lambda is such that the mean is matched *)
		If[Length[initial]==0,
			alphav = Table[1./M,{M}];
			lambda = DiagonalMatrix[orders].Range[M];
			trm = Mean[trace];
			inim = Total[alphav/Range[M]];
			lambda = lambda * inim / trm;
			P = Table[1,{M},{1}].{alphav};
		, If[Length[initial]==2,
			lambda = initial[[1]];
			P = initial[[2]];
			alphav = DTMCSolve[[P]];
		,
			Throw["MAPFromTrace: Invalid initial branch probability and rate vectors!"];
		]];
    
		Ascale = Table[0.,{K}];
		Bscale = Table[0.,{K}];
		ologli = 1;
		logli = 0.00000000001;
		steps = 1;
		t1 = AbsoluteTime[];
		While[Abs[(ologli-logli)/logli]>stopCond && steps<=maxIter,
			ologli = logli;
			(* E-step: *)
			(* branch densities: *)
			Q={};
			Do[AppendTo[Q,((lambda[[i]]*trace)^(orders[[i]]-1) / (orders[[i]]-1)! * lambda[[i]]) * Exp[-lambda[[i]]*trace]],{i,M}];
			(* forward likelihood vectors: *)
			prev = alphav;
			scprev = 0;
			A={};
			Do[
				prev = prev.DiagonalMatrix[Q[[All,k]]].P;
				scale = Log2[Total[prev]];
				prev = prev * 2^-scale;
				Ascale[[k]] = scprev + scale;
				AppendTo[A, prev];
				scprev = Ascale[[k]];
			,{k,K}];
			Av = Prepend[A[[1;;-2,All]],alphav];
			Ascalev = Prepend[Ascale[[1;;-2]],0];
			(* backward likelihood vectors: *)
			next = Table[1,{M}];
			scprev = 0;
			B = {};
			Do[
				next = DiagonalMatrix[Q[[All,k]]].P.next;
				scale = Log2[Total[next]];
				next = next * 2^-scale;
				Bscale[[k]] = scprev + scale;
				PrependTo[B,next];
				scprev = Bscale[[k]];
			,{k,K,1,-1}];
			Bv = Append[B[[2;;,All]],Table[1,{M}]];
			Bscalev = Append[Bscale[[2;;]], 0];
 
			llh = alphav.B[[1,All]];
			logli = (Log[llh] + Bscale[[1]] * Log[2]) / K;
			illh = 1.0 / llh;

			(* M-step: *)
			(* Calculate new estimates for the parameters *)
			AB = Av*B;
			nor = Total[AB,{2}];
			Do[AB[[All,m]] = AB[[All,m]]/nor, {m,M}];
			v1 = Total[AB,{1}];
			v2 = trace.AB;
			alphav = v1/K;
			lambda = orders*v1/v2;
        
			Avv = Av*Transpose[Q];
			nor = illh*2^(Ascalev+Bscalev-Bscale[[1]]);
			Do[Avv[[All,m]]=Avv[[All,m]]*nor,{m,M}];   
			P = (Transpose[Avv].Bv)*P;
			Do[P[[m,All]] = P[[m,All]] / Total[P[[m,All]]], {m,M}];

			steps++;
			If[BuTools`Verbose && AbsoluteTime[]-t1>2,
				Print["Num of iterations: ",steps, ", logli: ",logli];
				t1 = AbsoluteTime[];
			];
		];

		If[BuTools`Verbose,
			Print["Num of iterations: ",steps,", logli: ", logli];
			Print["EM algorithm terminated. (orders=",orders,")"];
        ];

		If[result=="vecmat",
			D0 = lambda;
			D1 = P;
		, If[result=="matmat",
			NN = Total[orders];
			D0 = Table[0,{NN},{NN}];
			ix = 1;
			Do[
				blk = -DiagonalMatrix[Table[1,{orders[[i]]}]];
				If[orders[[i]]>1, blk+=DiagonalMatrix[Table[1,{orders[[i]]-1}],1]];
				D0[[ix;;ix+orders[[i]]-1, ix;;ix+orders[[i]]-1]] = lambda[[i]]*blk;
				ix += orders[[i]];
			,{i,M}];
			D1 = Table[0.,{NN},{NN}];
			indicesTo = Prepend[Accumulate[orders[[1;;-2]]]+1,1];
			indicesFrom = Accumulate[orders];
			D1[[indicesFrom,indicesTo]] = DiagonalMatrix[lambda].P;
		]];
	];
    Return[{D0,D1,logli}];
];


End[(* Private *)];
EndPackage[];
