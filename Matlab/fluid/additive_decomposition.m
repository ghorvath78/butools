function [pmatrix, f, x]=additive_decomposition(Q,R,thres)
%Solves the finite fluid queue if the system is homogeneous between
%threshold levels
%based on the article: Nail Akar - Solving Multi-Regime fluid Queues
%it assumes that the rate is nonzero for all elements of the states
%space

%---------------------------notes------------------------------------------
%thres is a vector which contains the threshold levels of the system. Q is
%the generator matrix. Q(:,:,i) is the generator matrix for the i-th
%regime. Similarily R(:,:,i) is a diagonal matrix containg the rates for
%the i-th regime.

%the probabilities at the threshold levels are in the variable pmatrix
%the values of the density for the different states will be stored in the
%variable f
n=100;    %the number of points for which it calculates the probability density
         %as the program uses analytical normalization a low n will not result in
         %inaccuracy

%---------------------------parameters-------------------------------------
%thres is a row vector which contains the threshold levels. It's first
%element is the lower buffer limit (usually 0), the last element is the upper
%buffer limit (denoted with B in the article)
%if thres is a coloumn-vector we transform it to a rowvector
if size(thres,2)==1;
    thres=thres';
end
thresnum=size(thres,2); %the number of thresholds (containing the upper and lower buffer limit)

statenum=size(Q,1); %the cardinality of the state space of the Markov-chain
%Q contains the generatormatrix for every homogeneous part. It's size is (statenum,statenum,thresnum-1)
%R is the diagonal matrix with the size (statenum,statenum,thresnum-1)
%containg the (nonzero) rates for every state in every homogeneous part

A=zeros(size(Q)); %By multiplying both sides of the differential equation with R^(-1) we obtain A. It is LOWER level-dependent
for l=1:thresnum-1
    A(:,:,l)=Q(:,:,l)/R(:,:,l);
end
%-------------------------The additive decomposition-----------------------
%The matrix A is decomposed using the generalized Schur form
disp('Additive decomposition')
%The transformational matrices  of equation (23)
Y=zeros(size(Q)); 
Yinv=Y;
sumq=zeros(thresnum-1,1);
sump=sumq;

for l=1:thresnum-1 %a decomposition is needed on all levels
    tempA=A(:,:,l);
    [TRANSZ, schurA]=schur(tempA);
    
    %rearraging the blocks based on the sign eigenvalues (first the 0, than the positive eigenvalues, followed by the negative ones)
    q=0; p=0; %the number of negative and positive eigenvalues
    a=ordereigs(schurA); %for new Matlabs ordeig() is an availble function, which can replace this function (just replace ordereigs with ordeig)
    zeroeig=find(abs(a)==min(abs(a))); %the coordinate of the state with the zero eigenvalue
    
    pos=zeros(size(a)); %the coordinates positive eigenvalues will be stores in this vector
    neg=pos;
    temp=pos; %this will only contain the coordinate of the zero eigenvalue
    temp(zeroeig)=3;
    
    for j=1:statenum
        if (a(j)>0 && j~=zeroeig)
            pos(j)=1;
            p=p+1;
        elseif (a(j)<0 && j~=zeroeig)
           neg(j)=2;
           q=q+1;
        end
    end

    sor=pos+neg+temp; %this is a vector whose i-th element is 1, if schurA(i,i)>0, 2 if <0 and 3 if ==0
    
    %in matlab ordschur is a built-in function,  but it is unavailable for
    %octave
    [TRANSZ, schurA]=ordschur(TRANSZ,schurA,sor);  %the new schurA will be ordered according to sor    	       	
    
    %additive decomposition using appendix B
    X1=-schurA(1,2:end)/schurA(2:end,2:end);
    X2=lyapunov(-schurA(2:1+q,2:1+q),schurA(2+q:end,2+q:end),schurA(2:1+q,2+q:end));

    I=eye(statenum,statenum);
    transz1=I;
    transz1(1:size(X1,1),end-size(X1,2)+1:end)=-X1;
    transz2=I;
    transz2(end-size(X2,2)+1-size(X2,1):end-size(X2,2),end-size(X2,2)+1:end)=-X2;
    
    Y(:,:,l)=TRANSZ*transz1*transz2;
    Yinv(:,:,l)=Y(:,:,l)\I;
    sumq(l)=q;
    sump(l)=p;
   
end

decomposedA=zeros(size(Q));
for l=1:thresnum-1
decomposedA(:,:,l)=Yinv(:,:,l)*A(:,:,l)*Y(:,:,l);


%checking, if the decomposition worked
A1=decomposedA(1,1,l);
A2=decomposedA(2:sumq(l)+1,2:sumq(l)+1,l);
A3=decomposedA(sumq(l)+2:end,sumq(l)+2:end,l);
TEMP=[A1, zeros(1,size(A2,2)+size(A3,2)); 
      zeros(size(A2,1),1), A2, zeros(size(A2,1),size(A3,2));
      zeros(size(A3,1),1+size(A2,2)),A3];
  
    if norm(TEMP-decomposedA(:,:,l))>norm(TEMP)*1e-6
        warning('unable to decompose the matrix for all intervals')
    end
end

%---------------------------Linear equations-------------------------------
disp('creating the linear equations')
%Creating the matrix containg the boundary conditions
%The solution will be obtained by obtaining the rowvector pmatrix which
%containes the provababilities at the threshold levels. The a0,a-,a+
%values at the threshold levels will be stores in the prds vector

%the equations will take the form: pmatrix*A1+prds*A2=0
A1=[];
A2=[];

%the equations at the lower buffer-limit
Yinvtemp=Yinv(:,:,1);
Y0=Yinvtemp(1,:);
Ym=Yinvtemp(2:1+sumq(1),:);
Yp=Yinvtemp(2+sumq(1):end,:);
Am=decomposedA(2:sumq(1)+1,2:sumq(1)+1,1); %#ok<*NASGU>
Ap=decomposedA(2+sumq(1):end,2+sumq(1):end,1);

TEMP=[Y0;Ym;expm(-Ap*(thres(2)-thres(1)))*Yp]*R(:,:,1);

for l=1:statenum
    equation1=zeros(statenum*thresnum,1);
    equation2=zeros(statenum*(thresnum-1),1);
    if R(l,l,1)>0 %states with positive rate -> the probability is 0 at this level
        equation1(l)=1;
        A1=[A1,equation1]; %#ok<AGROW>
        %equation2 remains 0
        A2=[A2,equation2]; %#ok<AGROW>
    end
    
    %the second equation
    equation1(1:statenum)=-Q(:,l,1);
    equation2(1:statenum)=TEMP(:,l);
   
    A1=[A1,equation1]; %#ok<AGROW>
    A2=[A2,equation2]; %#ok<AGROW>
end

%the equations at the other threshold levels except the upper buffer limit
for m=2:thresnum-1 %for every threshold level
    Yinvtemp=Yinv(:,:,m-1);%the Yinv just below the threshold
    Y01=Yinvtemp(1,:);
    Ym1=Yinvtemp(2:1+sumq(m-1),:);
    Yp1=Yinvtemp(2+sumq(m-1):end,:);
    Am=decomposedA(2:sumq(m-1)+1,2:sumq(m-1)+1,m-1);
    
    Yinvtemp=Yinv(:,:,m);%the Yinv just above the threshold level
    Y02=Yinvtemp(1,:);
    Ym2=Yinvtemp(2:1+sumq(m),:);
    Yp2=Yinvtemp(2+sumq(m):end,:);
    Ap=decomposedA(2+sumq(m):end,2+sumq(m):end,m);
    
    LOWER=[Y01;expm(Am*(thres(m)-thres(m-1)))*Ym1;Yp1];
    UPPER=[Y02; Ym2; expm(-Ap*(thres(m+1)-thres(m)))*Yp2];
     
    LOWER=LOWER*R(:,:,m-1);
    LOWERTEMP=LOWER;
    UPPER=UPPER*R(:,:,m);
      
    for l=1:statenum
        equation1=zeros(statenum*thresnum,1);
        equation2=zeros(statenum*(thresnum-1),1);
        if R(l,l,m-1)>0 && R(l,l,m)>0
            %the probability at the threshold level is 0
            equation1((m-1)*statenum+l)=1;
            A1=[A1,equation1]; %#ok<AGROW>
            %equation2 remains 0
            A2=[A2,equation2]; %#ok<AGROW
       %elseif R(l,l,m-1)>0 && R(l,l,m)<0
            %there is no equation in this condition
        elseif R(l,l,m-1)<0 && R(l,l,m)>0
            %the probability at the threshold level is 0
            equation1((m-1)*statenum+l)=1;
            A1=[A1,equation1]; %#ok<AGROW>
            %equation2 remains 0
            A2=[A2,equation2]; %#ok<AGROW>
            %the density function is 0 just below the threshold level
            equation1((m-1)*statenum+l)=0;%this step is necessary because this was set to 1 for thr previous equation
            A1=[A1,equation1]; %#ok<AGROW>
            equation2((m-2)*statenum+1:(m-1)*statenum)=LOWERTEMP(:,l);
            A2=[A2, equation2]; %#ok<AGROW>
        elseif R(l,l,m-1)<0 && R(l,l,m)<0
            %the probability at the threshold level is 0
            equation1((m-1)*statenum+l)=1;
            A1=[A1,equation1]; %#ok<AGROW>
            %equation2 remains 0
            A2=[A2,equation2]; %#ok<AGROW>
        end
        %the following equation is true for every state
        equation1=zeros(statenum*thresnum,1);
        equation2=zeros(statenum*(thresnum-1),1);
        
        equation1((m-1)*statenum+1:m*statenum)=Q(:,l,m); %right continuous
        A1=[A1,equation1]; %#ok<AGROW>
        
        equation2((m-2)*statenum+1:(m-1)*statenum)=LOWER(:,l);
        equation2((m-1)*statenum+1:m*statenum)=-UPPER(:,l);
        A2=[A2,equation2]; %#ok<AGROW>
    end
end

%Linear equations at the upper buffer limit
Yinvtemp=Yinv(:,:,end);
Y0=Yinvtemp(1,:);
Ym=Yinvtemp(2:1+sumq(end),:);
Yp=Yinvtemp(2+sumq(end):end,:);
Am=decomposedA(2:sumq(end)+1,2:sumq(end)+1,end);

TEMP=[Y0;expm(Am*(thres(end)-thres(end-1)))*Ym;Yp]*R(:,:,end);
    
for l=1:statenum
    equation1=zeros(statenum*thresnum,1);
    equation2=zeros(statenum*(thresnum-1),1);
    if R(l,l,end)<0 %the probability is 0 for states with negative drift
        equation1((thresnum-1)*statenum+l)=1;
        A1=[A1,equation1]; %#ok<AGROW>
        %equation2 does not change
        A2=[A2,equation2]; %#ok<AGROW>
    end

    %the second equation
    equation1((thresnum-1)*statenum+1:end)=Q(:,l,end);
    equation2((thresnum-2)*statenum+1:end)=TEMP(:,l);
    
    A1=[A1,equation1]; %#ok<AGROW>
    A2=[A2,equation2]; %#ok<AGROW>
end

%solving the linear equations
disp('solving...')
LINEQ=[A1;A2];
sol=null(LINEQ')';

hely=0;
if rank(LINEQ)~=size(LINEQ,1)-1
    warning('possible numerical instability')
    [sv,se]=eig(LINEQ');
    se=abs(diag(se));
    hely=find(se(:)==min(se(:)));
    sol=sv(:,hely)';
    sol=sol(1,:);
end

%-----------------------obtaining the data, normalization------------------
pmatrix=sol(1:thresnum*statenum)';
a=sol(thresnum*statenum+1:end)';

p=zeros(thresnum,statenum);
temp=zeros(thresnum-1,statenum); %matrix, with rows made up of the a-s
for l=1:thresnum
    p(l,:)=pmatrix((l-1)*statenum+1:l*statenum);
    if l<thresnum
        temp(l,:)=a((l-1)*statenum+1:l*statenum);
    end
end

%setting the signs of p:
if norm(p-abs(p))>1e-6
    p=-p;
    temp=-temp;
end

%checking, if all values of p are nonnegative, it usually occurs when there
%are two solutions to LINEQ
if norm(p-abs(p))>1e-6 
    if size(hely,1)>1
        sol=sv(:,hely(2))';  
        for l=1:thresnum
            p(l,:)=pmatrix((l-1)*statenum+1:l*statenum);
            if l<thresnum
                temp(l,:)=a((l-1)*statenum+1:l*statenum);
            end
        end
        %setting the signs of p:
        if norm(p-abs(p))>1e-6
            p=-p;
            temp=-temp;
        end
    end
    if norm(p-abs(p))>1e-6
            warning('some elements of p are negative');
    end
end


%normalization
possz=sum(sum(p)); %the sum of the probability vector
fossz=0; %the sum of the probability density
for l=2:thresnum
    Yinvtemp=Yinv(:,:,l-1);
    Y0=Yinvtemp(1,:);
    Ym=Yinvtemp(2:1+sumq(l-1),:);
    Yp=Yinvtemp(2+sumq(l-1):end,:);
    Am=decomposedA(2:sumq(l-1)+1,2:sumq(l-1)+1,l-1);
    Ap=decomposedA(2+sumq(l-1):end,2+sumq(l-1):end,l-1);
    Aminv=Am\eye(sumq(l-1)); Apinv=Ap\eye(sump(l-1));
    TEMP=[Y0.*(thres(l)-thres(l-1));
          Aminv*(expm(Am.*(thres(l)-thres(l-1)))-eye(size(Am)))*Ym;
          Apinv*(eye(size(Ap))-expm(-Ap.*(thres(l)-thres(l-1))))*Yp];
    fossz=fossz+sum(temp(l-1,:)*TEMP);
end
normalizator=fossz+possz;

sfgv=temp/normalizator;
pmatrix=p/normalizator;


%----------------------------post processing-------------------------------
disp('accuracy:')
disp(min(norm(sol*LINEQ)))
disp('probabilities at threshold levels:')
disp('   level    probability')
disp([thres(:),sum(pmatrix,2)])

h=(thres(end)-thres(1))/(n-1);
k=0;
f=zeros(n,statenum);
for x=0:h:thres(end)
    k=k+1;
    for l=1:thresnum-1
        if x==thres(l)
            f(k,:)=NaN;
        elseif thres(l)<x && x<thres(l+1)
            Yinvtemp=Yinv(:,:,l);
            Y0=Yinvtemp(1,:);
            Ym=Yinvtemp(2:1+sumq(l),:);
            Yp=Yinvtemp(2+sumq(l):end,:);
            Am=decomposedA(2:sumq(l)+1,2:sumq(l)+1,l);
            Ap=decomposedA(2+sumq(l):end,2+sumq(l):end,l);
            TEMP=[Y0;
                  expm(Am.*(x-thres(l)))*Ym;
                  expm(-Ap.*(thres(l+1)-x))*Yp];
             f(k,:)=sfgv(l,:)*TEMP;
        end
    end 
end
f(end,:)=NaN;

x=0:h:thres(end);
end








