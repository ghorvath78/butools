function [PSZI, PSZI2, lambdapp, lambdamm, pszipm, pszimp] = matrix_analytic_psi(E, T, puff)
%It calculates PSZI and other similar matrices

statenum=size(T,1);
%first it orders the states corresponding to their rates (0 rate, positive
%rate then negative rate)
ORDERMATRIX=zeros(statenum,statenum);
nullnum=0; poz=0; neg=0;
zerovector=zeros(statenum,1); posvector=zeros(statenum,1); negvector=zeros(statenum,1);
for i=1:statenum;
    if E(i,i)==0
       nullnum=nullnum+1; %the number of state with 0 rate
       zerovector(nullnum)=i;       
    elseif E(i,i)>0;
        poz=poz+1; %the number of states with>0 rate
        posvector(poz)=i;
    else
        neg=neg+1;
        negvector(neg)=i;
    end
end
        
for i=1:statenum;
    if i<=nullnum;
        ORDERMATRIX(zerovector(i),i)=1;
    elseif i<=nullnum+poz
        ORDERMATRIX(posvector(i-nullnum),i)=1;
    else
        ORDERMATRIX(negvector(i-nullnum-poz),i)=1;
    end
end


E=ORDERMATRIX'*E*ORDERMATRIX;
T=ORDERMATRIX'*T*ORDERMATRIX;

T=abs(E)^(-1)*T;


%-------------------------calculating G------------------------------------
I=eye(statenum);
Ap=zeros(statenum);
A0=zeros(statenum);
Am=zeros(statenum);

c=max(abs(diag(T))); 

P=I+1/c*T;

if (poz+neg)~=statenum
    warning('this matrix should not have states with 0 rate');
end

Ap(1:poz,1:poz)=eye(poz,poz)/2;
A0(1:poz,1:poz)=P(1:poz,1:poz)/2;
A0(poz+1:poz+neg, 1:poz)=P(poz+1:poz+neg,1:poz);
Am(1:poz,poz+1:poz+neg)=1/2*P(1:poz,poz+1:poz+neg);
Am(poz+1:poz+neg,poz+1:poz+neg)=P(poz+1:poz+neg,poz+1:poz+neg);

A2=Ap;
A1=A0;
A0=Am;

if max(sum(A0+A1+A2,2)-ones(size(A0,1),1))>1e-6
    warning('A0+A1+A2 is not (sub)stochastic')
end

G=QBD_CR(A0,A1,A2);

PSZI=G(1:poz, poz+1:neg+poz);

%--------------------------------------------------------------------------
%the same procedure for the reversed process (we switch the sign of the
%rates)

%we need to transform the matrices so that we can use the same algorithm as
%before
T2(1:neg,1:neg)=T(poz+1:poz+neg,poz+1:poz+neg);
T2(1:neg,neg+1:neg+poz)=T(poz+1:poz+neg,1:poz);
T2(neg+1:poz+neg,neg+1:poz+neg)=T(1:poz,1:poz);
T2(neg+1:poz+neg,1:neg)=T(1:poz,poz+1:poz+neg);


I=eye(size(T));
Ap=zeros(size(T));
A0=zeros(size(T));
Am=zeros(size(T));

c=max(abs(diag(T2))); 

P=I+1/c*T2;

%the number of states with positive and negative rate
temp=neg;
neg=poz;
poz=temp;


Ap(1:poz,1:poz)=eye(poz,poz)/2;
A0(1:poz,1:poz)=P(1:poz,1:poz)/2;
A0(poz+1:poz+neg, 1:poz)=P(poz+1:poz+neg,1:poz);
Am(1:poz,poz+1:poz+neg)=1/2*P(1:poz,poz+1:poz+neg);
Am(poz+1:poz+neg,poz+1:poz+neg)=P(poz+1:poz+neg,poz+1:poz+neg);

A2=Ap;
A1=A0;
A0=Am;


if max(sum(A0+A1+A2,2))>1+1e-6
    disp('A0+A1+A2 is not (sub)stochastic')
end

G=QBD_CR(A0,A1,A2);

PSZI2=G(1:poz, poz+1:neg+poz);

%we set back the values of poz and the neg, to their original values
temp=neg;
neg=poz;
poz=temp;

%--------------------------------------------------------------------------
%calculation of the Lambda matrices

%solving the linear equationss>

U=T(poz+1:poz+neg,poz+1:poz+neg)+T(poz+1:poz+neg, 1:poz)*PSZI;
U2=T(1:poz,1:poz)+T(1:poz,poz+1:poz+neg)*PSZI2;
K=T(1:poz,1:poz)+PSZI*T(poz+1:poz+neg, 1:poz); %#ok<NASGU>
K2=T(poz+1:poz+neg,poz+1:poz+neg)+PSZI2*T(1:poz,poz+1:poz+neg); %#ok<NASGU>


lambdapp=((eye(size(PSZI*PSZI2))-PSZI*PSZI2)*expm(U2*puff)*(eye(size(PSZI*expm(U*puff)*PSZI2*expm(U2*puff)))-PSZI*...
    expm(U*puff)*PSZI2*expm(U2*puff))^(-1));
lambdamm=((eye(size(PSZI2*PSZI))-PSZI2*PSZI)*expm(U*puff)*(eye(size(PSZI2*expm(U2*puff)*PSZI*expm(U*puff)))-PSZI2*...
    expm(U2*puff)*PSZI*expm(U*puff))^(-1));


pszipm=(PSZI-expm(U2*puff)*PSZI*expm(U*puff))*(eye(size(PSZI2*expm(U2*puff)*PSZI*expm(U*puff)))-PSZI2*expm(U2*puff)*PSZI*expm(U*puff))^(-1);
pszimp=(PSZI2-expm(U*puff)*PSZI2*expm(U2*puff))*(eye(size(PSZI*expm(U*puff)*PSZI2*expm(U2*puff)))-PSZI*expm(U*puff)*PSZI2*expm(U2*puff))^(-1);