function [pmatrix, f, x]=matrix_analytic(Q,R,thres)
%Solves the finite fluid queue if the system is homogeneous between
%threshold levels
%based on the article: Ana da Silva Soares, Guy Latouche - Fluid queues with level dependent evolution
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
n=100;    %the number of points -1 for which it calculates the probability density
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

%The function can not handle systems with thresnum==2 or 3. Therefore we
%add two more thresholds. (This will not change the results, but may slow
%down the time it takes to solve the problem.)
Rorig=R;
thresnumtemp=0;
if thresnum==2 || thresnum==3
    thres=[thres(1),thres,thres(end)];
    thresnumtemp=thresnum;
    thresnum=thresnum+2;
    Qtemp=Q;
    Rtemp=R;
    Q=zeros(statenum,statenum,thresnum-1);
    R=Q;
    Q(:,:,1)=Qtemp(:,:,1);
    Q(:,:,2:end-1)=Qtemp(:,:,1:end);
    Q(:,:,end)=Qtemp(:,:,end);
    
    R(:,:,1)=Rtemp(:,:,1);
    R(:,:,2:end-1)=Rtemp(:,:,1:end);
    R(:,:,end)=Rtemp(:,:,end);
end


%----------------------calculation of PSZI----------------------------------
disp('calculation of PSZI')
%these will contain the corresponding matrices
PSZI=zeros(statenum,statenum,thresnum-1);
PSZI2=PSZI; lambdapp=PSZI; lambdamm=PSZI; pszipm=PSZI; pszimp=PSZI;

%the sizes of the matrices will be stored in these matrices
PSZIsize=zeros(thresnum-1,2); PSZI2size=zeros(thresnum-1,2);
lambdappsize=zeros(thresnum-1,2); lambdammsize=zeros(thresnum-1,2);
pszipmsize=zeros(thresnum-1,2); pszimpsize=zeros(thresnum-1,2);


for l=1:thresnum-1
    [PSZItemp, PSZI2temp, lambdapptemp, lambdammtemp, pszipmtemp, pszimptemp]=matrix_analytic_psi(R(:,:,l), Q(:,:,l), thres(l+1)-thres(l));
    PSZIsize(l,:)=size(PSZItemp);
    PSZI(1:PSZIsize(l,1),1:PSZIsize(l,2),l)=PSZItemp;
    PSZI2size(l,:)=size(PSZI2temp);
    PSZI2(1:PSZI2size(l,1),1:PSZI2size(l,2),l)=PSZI2temp;
    lambdappsize(l,:)=size(lambdapptemp);
    lambdapp(1:lambdappsize(l,1),1:lambdappsize(l,2),l)=lambdapptemp;
    lambdammsize(l,:)=size(lambdammtemp);
    lambdamm(1:lambdammsize(l,1),1:lambdammsize(l,2),l)=lambdammtemp;
    pszipmsize(l,:)=size(pszipmtemp);
    pszipm(1:pszipmsize(l,1),1:pszipmsize(l,2),l)=pszipmtemp;
    pszimpsize(l,:)=size(pszimptemp);
    pszimp(1:pszimpsize(l,1),1:pszimpsize(l,2),l)=pszimptemp;
    
end

for l=1:thresnum-1
    PSZI(:,:,l)=sparse(PSZI(:,:,l));
    PSZI2(:,:,l)=sparse(PSZI2(:,:,l));
    lambdapp(:,:,l)=sparse(lambdapp(:,:,l));
    lambdamm(:,:,l)=sparse(lambdamm(:,:,l));
    pszipm(:,:,l)=sparse(pszipm(:,:,l));
    pszimp(:,:,l)=sparse(pszimp(:,:,l));
end
  

%-------------------calculation of the probabilities-----------------------
disp('calculation of the probabilities')

TRANSmatrix=zeros(statenum,statenum,thresnum);
upnumvector=zeros(thresnum,1);
sticknumvector=zeros(thresnum,1);
downnumvector=zeros(thresnum,1);
repulsivenumvector=zeros(thresnum,1);

%--------------------mátrices a the lower buffer limit---------------------

%-----------ordering based on their types (sticky/repulsive etc.)----------
I=eye(statenum,statenum);
%there are only sticky and repulsive states at the bottom

P=I+diag(diag(-Q(:,:,1)))\I*Q(:,:,1);
%the states are separated based on their corresponding rates
pos=0;  posvector=zeros(statenum,1);
neg=0;  negvector=zeros(statenum,1);
TRANS=zeros(statenum,statenum);
for l=1:statenum
    if R(l,l,1)<0
        neg=neg+1;
        negvector(neg)=l;
    else
        pos=pos+1;
        posvector(pos)=l;
    end
end

for l=1:statenum
    if l<=pos
        TRANS(posvector(l),l)=1;
    else
        TRANS(negvector(l-pos),l)=1;
    end
end

%Saving the matrices, as they will be used later
TRANSmatrix(:,:,1)=TRANS;
upnumvector(1)=pos;
sticknumvector(1)=neg;
downnumvector(1)=0;
repulsivenumvector(1)=0;


P=TRANS'*P*TRANS; P=P(pos+1:end,:);
Pszius0=pszipm(1:pszipmsize(1,1),1:pszipmsize(1,2),1);
Pss0=P(:,pos+1:end);
Psu0=P(:,1:pos);

%calculation of lambadauu0, lambdaus0
indexvector=[];
for l=1:statenum
    if R(l,l,1)>0 && R(l,l,2)>0 %the first threshold is up
        indexvector=[indexvector;1];
    elseif R(l,l,1)>0 && R(l,l,2)<0 %sticky
        indexvector=[indexvector;2];
    end
end

upnum=0; sticknum=0;
upvector=zeros(statenum,1);
stickvector=zeros(statenum,1);
for l=1:size(indexvector,1)
    if indexvector(l)==1
        upnum=upnum+1;
        upvector(upnum)=l;
    else% if indexvector==2
        sticknum=sticknum+1;
        stickvector(sticknum)=l;
    end
end

TRANS=zeros(size(indexvector,1));
for l=1:size(indexvector,1)
    if l<= upnum
        TRANS(upvector(l),l)=1;
    else
        TRANS(stickvector(l-upnum),l)=1;
    end
end

lambdapptemp=lambdapp(1:lambdappsize(1,1),1:lambdappsize(1,2),1);
lambdapptemp=lambdapptemp*TRANS;

lambdauu0=lambdapptemp(:,1:upnum);
lambdaus0=lambdapptemp(:,upnum+1:end);
        
%pszidu0, pszids0
pszi2temp=pszimp(1:pszimpsize(1,1),1:pszimpsize(1,2),1);
%the ordering of the coloumn is the same:
pszi2temp=pszi2temp*TRANS;

%when ordering the rows, we should remove the repulsive states at thres(2)

indexvector=[];
for l=1:statenum
    if R(l,l,1)<0 && R(l,l,2)<0 %down
        indexvector=[indexvector;3];
    elseif R(l,l,1)<0 && R(l,l,2)>0 %repulsive
        indexvector=[indexvector;4];
    end
end

downnum=0; repulsivenum=0;
downvector=zeros(statenum,1);
repulsivevector=zeros(statenum,1);
for l=1:size(indexvector,1)
    if indexvector(l)==3
        downnum=downnum+1;
        downvector(downnum)=l;
    else%if indexvector==4
        repulsivenum=repulsivenum+1;
        repulsivevector(repulsivenum)=l;
    end
end

TRANS2=zeros(size(indexvector,1));
for l=1:size(indexvector,1)
    if l<= downnum
        TRANS2(downvector(l),l)=1;
    else
        TRANS2(repulsivevector(l-downnum),l)=1;
    end
end

pszi2temp=TRANS2'*pszi2temp;

pszidu0=pszi2temp(1:downnum,1:upnum);
pszids0=pszi2temp(1:downnum,upnum+1:end);

%calculation of lambdads0
%we will need to remove the repulsive states with the use of TRANS2
lambdads0=lambdamm(1:lambdammsize(1,1),1:lambdammsize(1,2),1);
lambdads0=TRANS2'*lambdads0;
lambdads0=lambdads0(1:downnum,:);

%----------------------------saving----------------------------------------
psziussize=zeros(thresnum,2);
psziudsize=zeros(thresnum,2);
pszirssize=zeros(thresnum,2);
pszirdsize=zeros(thresnum,2);
psusize=zeros(thresnum,2);
psrsize=zeros(thresnum,2);
psssize=zeros(thresnum,2);
psdsize=zeros(thresnum,2);
lambdauusize=zeros(thresnum,2);
lambdaussize=zeros(thresnum,2);
lambdarusize=zeros(thresnum,2);
lambdarssize=zeros(thresnum,2);
lambdadssize=zeros(thresnum,2);
lambdaddsize=zeros(thresnum,2);
pszidusize=zeros(thresnum,2);
pszidssize=zeros(thresnum,2);

pszius=zeros(statenum,statenum,thresnum);
psziud=zeros(statenum,statenum,thresnum);
pszirs=zeros(statenum,statenum,thresnum);
pszird=zeros(statenum,statenum,thresnum);
psu=zeros(statenum,statenum,thresnum);
psr=zeros(statenum,statenum,thresnum);
pss=zeros(statenum,statenum,thresnum);
psd=zeros(statenum,statenum,thresnum);
lambdauu=zeros(statenum,statenum,thresnum);
lambdaus=zeros(statenum,statenum,thresnum);
lambdaru=zeros(statenum,statenum,thresnum);
lambdars=zeros(statenum,statenum,thresnum);
lambdads=zeros(statenum,statenum,thresnum);
lambdadd=zeros(statenum,statenum,thresnum);
pszidu=zeros(statenum,statenum,thresnum);
pszids=zeros(statenum,statenum,thresnum);

%saving the data
psziussize(1,:)=size(Pszius0);
pszius(1:psziussize(1,1),1:psziussize(1,2),1)=Pszius0;
psssize(1,:)=size(Pss0);
pss(1:psssize(1,1),1:psssize(1,2),1)=Pss0;
psusize(1,:)=size(Psu0);
psu(1:psusize(1,1),1:psusize(1,2),1)=Psu0;
lambdadssize(1,:)=size(lambdads0);
lambdads(1:lambdadssize(1,1),1:lambdadssize(1,2),1)=lambdads0;
lambdauusize(1,:)=size(lambdauu0);
lambdauu(1:lambdauusize(1,1),1:lambdauusize(1,2),1)=lambdauu0;
lambdaussize(1,:)=size(lambdaus0);
lambdaus(1:lambdaussize(1,1),1:lambdaussize(1,2),1)=lambdaus0;
pszidusize(1,:)=size(pszidu0);
pszidu(1:pszidusize(1,1),1:pszidusize(1,2),1)=pszidu0;
pszidssize(1,:)=size(pszids0);
pszids(1:pszidssize(1,1),1:pszidssize(1,2),1)=pszids0;

%-------------------matrices at the intermediate levels--------------------
for l=2:thresnum-1
    %separating the sticky, repulsive etc. states for the lth threshold
    indexvector=zeros(statenum,1);
    for m=1:statenum
        if R(m,m,l-1)>0 && R(m,m,l)>0 %up
            indexvector(m)=1;
        elseif R(m,m,l-1)>0 && R(m,m,l)<0 %sticky
            indexvector(m)=2;
        elseif R(m,m,l-1)<0 && R(m,m,l)<0 %down
            indexvector(m)=3;
        else %repulsive
            indexvector(m)=4;
        end
    end
    
    %creating the ordering matrix:
     upnum=0;   upvector=zeros(statenum,1);
     sticknum=0;  stickvector=zeros(statenum,1);
     downnum=0; downvector=zeros(statenum,1);
     repulsivenum=0; repulsivevector=zeros(statenum,1);
     for m=1:statenum
        if indexvector(m)==1
            upnum=upnum+1;
            upvector(upnum)=m;
        elseif indexvector(m)==2
            sticknum=sticknum+1;
            stickvector(sticknum)=m;
        elseif indexvector(m)==3
            downnum=downnum+1;
            downvector(downnum)=m;
        else
            repulsivenum=repulsivenum+1;
            repulsivevector(repulsivenum)=m;
        end
     end
     TRANS=zeros(statenum);
     
     for m=1:statenum
         if m<=upnum
             TRANS(upvector(m),m)=1;
         elseif m<=upnum+sticknum
             TRANS(stickvector(m-upnum),m)=1;
         elseif m<=upnum+sticknum+downnum
             TRANS(downvector(m-upnum-sticknum),m)=1;
         else
             TRANS(repulsivevector(m-upnum-sticknum-downnum),m)=1;
         end
     end
     %saving:
     TRANSmatrix(:,:,l)=TRANS;
     upnumvector(l)=upnum;
     sticknumvector(l)=sticknum;
     downnumvector(l)=downnum;
     repulsivenumvector(l)=repulsivenum;
     
    %CREATING THE LAMBDAPP MATRICES   
    %the up and sticky states should also be separated for the l+1th 
    %threshold to calculate the lambdauu etc. matrices
    indexvectorupper=[];
    if l<thresnum-1
        for m=1:statenum
            if R(m,m,l)>0 && R(m,m,l+1)>0 %up
                indexvectorupper=[indexvectorupper;1];
            elseif R(m,m,l)>0 && R(m,m,l+1)<0 %sticky
                indexvectorupper=[indexvectorupper;2];
            end
        end
    else%for the highest thershold-level
        for m=1:statenum
            if R(m,m,l)>0%sticky
                indexvectorupper=[indexvectorupper;2];
            end
        end
    end
    
   upnumupper=0; upvector=zeros(statenum,1);
   sticknumupper=0; stickvector=zeros(statenum,1);
   for m=1:size(indexvectorupper,1)
        if indexvectorupper(m)==1
            upnumupper=upnumupper+1;
            upvector(upnumupper)=m;
        elseif indexvectorupper(m)==2
            sticknumupper=sticknumupper+1;
            stickvector(sticknumupper)=m;
        end
   end
   
   TRANS2=zeros(size(indexvectorupper,1));   
     for m=1:size(indexvectorupper,1)
         if m<=upnumupper
             TRANS2(upvector(m),m)=1;
         else
             TRANS2(stickvector(m-upnumupper),m)=1;
         end
     end
     
    
    %we need to order the matrices:
    lambdapptemp=lambdapp(1:lambdappsize(l,1),1:lambdappsize(l,2),l);
    %the ordering of coloumns
    lambdapptemp=lambdapptemp*TRANS2;
    %the ordering of rows
    TRANStemp=[TRANS(:,1:upnum), TRANS(:,upnum+sticknum+downnum+1:end)];
    theend=0; m=0;
     if size(TRANStemp,2)~=size(TRANStemp,1)
                while theend==0
                  m=m+1;
                  if norm(TRANStemp(m,:))<1e-10
                    TRANStemp=[TRANStemp(1:m-1,:);TRANStemp(m+1:end,:)];
                    m=m-1;
                    if size(TRANStemp,2)==size(TRANStemp,1)
                            theend=1;
                    end                
                  end
                end
     end
    
    %{
    %checking if the results are correct
    indexvectortemp=[];
    for m=1:statenum
        if indexvector(m)==1
            indexvectortemp=[indexvectortemp;1];
        elseif indexvector(m)==4
            indexvectortemp=[indexvectortemp;4];
        end
    end
    TRANStemp=zeros(size(indexvectortemp,1));
    oszlopindex=0;
    for m=1:size(indexvectortemp,1)
        if indexvectortemp(m)==1
            oszlopindex=oszlopindex+1;
            TRANStemp(m,oszlopindex)=1;
        end
    end

    for m=1:size(indexvectortemp,1)
        if indexvectortemp(m)==4
            oszlopindex=oszlopindex+1;
            TRANStemp(m,oszlopindex)=1;
        end
    end
    TRANStemp
    %}  
    
    lambdapptemp=TRANStemp'*lambdapptemp;
    
    %saving the data
    lambdauutemp=lambdapptemp(1:upnum,1:upnumupper);
    lambdarutemp=lambdapptemp(upnum+1:end,1:upnumupper);
    lambdaustemp=lambdapptemp(1:upnum,upnumupper+1:end);
    lambdarstemp=lambdapptemp(upnum+1:end,upnumupper+1:end);
    
    lambdauusize(l,:)=size(lambdauutemp);
    lambdauu(1:lambdauusize(l,1),1:lambdauusize(l,2),l)=lambdauutemp;
    lambdarusize(l,:)=size(lambdarutemp);
    lambdaru(1:lambdarusize(l,1),1:lambdarusize(l,2),l)=lambdarutemp;
    lambdaussize(l,:)=size(lambdaustemp);
    lambdaus(1:lambdaussize(l,1),1:lambdaussize(l,2),l)=lambdaustemp;
    lambdarssize(l,:)=size(lambdarstemp);
    lambdars(1:lambdarssize(l,1),1:lambdarssize(l,2),l)=lambdarstemp;
    
    
    %CREATING THE PSZI MATRICES
     pszitemp=pszipm(1:pszipmsize(l,1),1:pszipmsize(l,2),l);
     %ordering the coloumns based on the sticky and down elements of indexvector
     TRANStemp=[TRANS(:,upnum+1:upnum+sticknum), TRANS(:,upnum+sticknum+1:upnum+sticknum+downnum)];
     theend=0; m=0;
     if size(TRANStemp,2)~=size(TRANStemp,1)
                while theend==0
                  m=m+1;
                  if norm(TRANStemp(m,:))<1e-10
                    TRANStemp=[TRANStemp(1:m-1,:);TRANStemp(m+1:end,:)];
                    m=m-1;
                    if size(TRANStemp,2)==size(TRANStemp,1)
                            theend=1;
                    end                
                  end
                end
     end
     pszitemp=pszitemp*TRANStemp;
    
    %ordering the rows (up or repulsive)
    TRANStemp=[TRANS(:,1:upnum), TRANS(:,upnum+sticknum+downnum+1:end)];
    theend=0; m=0;
     if size(TRANStemp,2)~=size(TRANStemp,1)
                while theend==0
                  m=m+1;
                  if norm(TRANStemp(m,:))<1e-10
                    TRANStemp=[TRANStemp(1:m-1,:);TRANStemp(m+1:end,:)];
                    m=m-1;
                    if size(TRANStemp,2)==size(TRANStemp,1)
                            theend=1;
                    end                
                  end
                end
      end
     pszitemp=TRANStemp'*pszitemp;
     
     %saving the results
     psziustemp=pszitemp(1:upnum,1:sticknum);
     psziudtemp=pszitemp(1:upnum,sticknum+1:end);
     pszirstemp=pszitemp(upnum+1:end,1:sticknum);
     pszirdtemp=pszitemp(upnum+1:end,sticknum+1:end);
     
     psziussize(l,:)=size(psziustemp);
     pszius(1:psziussize(l,1),1:psziussize(l,2),l)=psziustemp;
     psziudsize(l,:)=size(psziudtemp);
     psziud(1:psziudsize(l,1),1:psziudsize(l,2),l)=psziudtemp;
     pszirssize(l,:)=size(pszirstemp);
     pszirs(1:pszirssize(l,1),1:pszirssize(l,2),l)=pszirstemp;
     pszirdsize(l,:)=size(pszirdtemp);
     pszird(1:pszirdsize(l,1),1:pszirdsize(l,2),l)=pszirdtemp;
     
    
     %CREATING THE P MATRICES
     P=I+diag(diag(-Q(:,:,l)))\I*Q(:,:,l);
     %all states (up,sticky, repulsive, down) are needed to be separated
     %using indexvector
     P=TRANS'*P*TRANS;
     P=P(upnum+1:upnum+sticknum,:);
     
     %saving the results
     psutemp=P(:,1:upnum);
     psstemp=P(:,upnum+1:upnum+sticknum);
     psdtemp=P(:,upnum+sticknum+1:upnum+sticknum+downnum);
     psrtemp=P(:,upnum+sticknum+downnum+1:end);
     
     psusize(l,:)=size(psutemp);
     psu(1:psusize(l,1),1:psusize(l,2),l)=psutemp;
     psssize(l,:)=size(psstemp);
     pss(1:psssize(l,1),1:psssize(l,2),l)=psstemp;
     psdsize(l,:)=size(psdtemp);
     psd(1:psdsize(l,1),1:psdsize(l,2),l)=psdtemp;
     psrsize(l,:)=size(psrtemp);
     psr(1:psrsize(l,1),1:psrsize(l,2),l)=psrtemp;
     
     %CREATING THE PSZI2 MATRICES
     %we also need to order the rows with TRANS2
     pszi2temp=pszimp(1:pszimpsize(l,1),1:pszimpsize(l,2),l);
     pszi2temp=pszi2temp*TRANS2;
    
     %the rows of pszi2temp are needed to be ordered, as they also contain 
     %repulsive states, not just down states
     %the upper indexvector is used for the separation of down and
     %repulsive states
     indexvectorupper2=[];
    if l<thresnum-1
        for m=1:statenum
            if R(m,m,l)<0 && R(m,m,l+1)<0 %down
                indexvectorupper2=[indexvectorupper2;3];
            elseif R(m,m,l)<0 && R(m,m,l+1)>0 %repulsive
                indexvectorupper2=[indexvectorupper2;4];
            end
        end
    else%the highest threshold
        for m=1:statenum
            if R(m,m,l)<0%down
                indexvectorupper2=[indexvectorupper2;3];
            end
        end
    end
    
   downnumupper=0; downvector=zeros(statenum,1);
   repulsivenumupper=0; repulsivevector=zeros(statenum,1);
   for m=1:size(indexvectorupper2,1)
        if indexvectorupper2(m)==3
            downnumupper=downnumupper+1;
            downvector(downnumupper)=m;
        elseif indexvectorupper2(m)==4
            repulsivenumupper=repulsivenumupper+1;
            repulsivevector(repulsivenumupper)=m;
        end
   end
   
   TRANS3=zeros(size(indexvectorupper2,1));   
     for m=1:size(indexvectorupper2,1)
         if m<=downnumupper
             TRANS3(downvector(m),m)=1;
         else
             TRANS3(repulsivevector(m-downnumupper),m)=1;
         end
     end
     
     pszi2temp=TRANS3'*pszi2temp;   
     
     %saving the data
     pszidutemp=pszi2temp(1:downnumupper,1:upnumupper);
     pszidstemp=pszi2temp(1:downnumupper,upnumupper+1:end);
     
     pszidusize(l,:)=size(pszidutemp);
     pszidu(1:pszidusize(l,1),1:pszidusize(l,2),l)=pszidutemp;
     pszidssize(l,:)=size(pszidstemp);
     pszids(1:pszidssize(l,1),1:pszidssize(l,2),l)=pszidstemp;
     
     %CREATING THE LAMBDAMM MATRICES
     lambdammtemp=lambdamm(1:lambdammsize(l,1),1:lambdammsize(l,2),l);
    
     %ordering the rows
     lambdammtemp=TRANS3'*lambdammtemp;
     %ordering the coloumns
     TRANStemp=[TRANS(:,upnum+1:upnum+sticknum), TRANS(:,upnum+sticknum+1:upnum+sticknum+downnum)];
     theend=0; m=0;
     if size(TRANStemp,2)~=size(TRANStemp,1)
                while theend==0
                  m=m+1;
                  if norm(TRANStemp(m,:))<1e-10
                    TRANStemp=[TRANStemp(1:m-1,:);TRANStemp(m+1:end,:)];
                    m=m-1;
                    if size(TRANStemp,2)==size(TRANStemp,1)
                            theend=1;
                    end                
                  end
                end
     end
     lambdammtemp=lambdammtemp*TRANStemp;

     %saving the results
     lambdadstemp=lambdammtemp(1:downnumupper,1:sticknum);
     lambdaddtemp=lambdammtemp(1:downnumupper,sticknum+1:end);
     
     lambdadssize(l,:)=size(lambdadstemp);
     lambdads(1:lambdadssize(l,1),1:lambdadssize(l,2),l)=lambdadstemp;
     lambdaddsize(l,:)=size(lambdaddtemp);
     lambdadd(1:lambdaddsize(l,1),1:lambdaddsize(l,2),l)=lambdaddtemp;
   
end

%-----------------------the matrices at the end----------------------------
%there are only two types of states at the end
P=I+diag(diag(-Q(:,:,end)))\I*Q(:,:,end); %here we use the fact, that Q does not change at the last trhreshold level


indexvector=zeros(statenum);
for l=1:statenum
    if R(l,l,end)>0; %stick
        indexvector(l)=2;
    else %down
        indexvector(l)=3;
    end
end

downnum=0; downvector=zeros(statenum,1);
sticknum=0; stickvector=zeros(statenum,1);
for m=1:statenum
      if indexvector(m)==2
            sticknum=sticknum+1;
            stickvector(sticknum)=m;
      elseif indexvector(m)==3
            downnum=downnum+1;
            downvector(downnum)=m;
      end
end

TRANS=zeros(statenum);   
for m=1:statenum
    if m<=sticknum
             TRANS(stickvector(m),m)=1;
    else
            TRANS(downvector(m-sticknum),m)=1;
    end
end
%saving the matrices
TRANSmatrix(:,:,end)=TRANS;
upnumvector(end)=0;
sticknumvector(end)=sticknum;
downnumvector(end)=downnum;
repulsivenumvector(end)=0;


P=TRANS'*P*TRANS;
P=P(1:sticknum,:);
psstemp=P(:,1:sticknum);
psdtemp=P(:,sticknum+1:end);

psssize(end,:)=size(psstemp);
pss(1:psssize(end,1),1:psssize(end,2),end)=psstemp;
psdsize(end,:)=size(psdtemp);
psd(1:psdsize(end,1),1:psdsize(end,2),end)=psdtemp;

%--------------------------------------------------------------------------
%-------------------constructing the omega matrix--------------------------
 %we create the omega matrix at the end of the article. S, U, L will be the
 %block matrices such that: omega=[S, U, ...;
 %                                 L, S, U ...;
 %                                ... L, S, U...]
 S=zeros(statenum,statenum,thresnum);
 U=zeros(statenum,statenum,thresnum-1);
 L=zeros(statenum,statenum,thresnum-1);
 
 %calculation S(:,:,1), U(:,:,1), L(:,:,1)

 S(:,:,1)=[zeros(psziussize(1,1),psusize(1,2)),pszius(1:psziussize(1,1),1:psziussize(1,2),1);
           psu(1:psusize(1,1),1:psusize(1,2),1),pss(1:psssize(1,1),1:psssize(1,2),1)];

 if lambdauusize(1,1) ~= 0
 U(1:lambdauusize(1,1),:,1)=[lambdauu(1:lambdauusize(1,1),1:lambdauusize(1,2),1),zeros(lambdauusize(1,1),psrsize(2,2)),lambdaus(1:lambdaussize(1,1),1:lambdaussize(1,2),1),zeros(lambdauusize(1,1),psdsize(2,2))];
else
 U(1:lambdaussize(1,1),:,1)=[lambdauu(1:lambdaussize(1,1),1:lambdauusize(1,2),1),zeros(lambdaussize(1,1),psrsize(2,2)),lambdaus(1:lambdaussize(1,1),1:lambdaussize(1,2),1),zeros(lambdaussize(1,1),psdsize(2,2))];  
end

L(end-lambdadssize(1,1)+1:end,end-lambdadssize(1,2)+1:end,1)=lambdads(1:lambdadssize(1,1),1:lambdadssize(1,2),1);

%calculating the matrices for the intermediate thresholds

for l=2:thresnum-1 
    
    
    S(:,:,l)=[zeros(max(psziussize(l,1),psziudsize(l,1)),psusize(l,2)+psrsize(l,2)),pszius(1:psziussize(l,1),1:psziussize(l,2),l),psziud(1:psziudsize(l,1),1:psziudsize(l,2),l);
              zeros(max(pszirssize(l,1),pszirdsize(l,1)),psusize(l,2)+psrsize(l,2)),pszirs(1:pszirssize(l,1),1:pszirssize(l,2),l),pszird(1:pszirdsize(l,1),1:pszirdsize(l,2),l);
              psu(1:psusize(l,1),1:psusize(l,2),l),psr(1:psrsize(l,1),1:psrsize(l,2),l),pss(1:psssize(l,1),1:psssize(l,2),l),psd(1:psdsize(l,1),1:psdsize(l,2),l);
              pszidu(1:pszidusize(l-1,1),1:pszidusize(l-1,2),l-1),zeros(max(pszidusize(l-1,1),pszidssize(l-1,1)),psrsize(l,2)),pszids(1:pszidssize(l-1,1),1:pszidssize(l-1,2),l-1),zeros(max(pszidusize(l-1,1),pszidssize(l-1,1)),psdsize(l,2))];
          
   
    U(1:lambdauusize(l,1)+lambdarusize(l,1),1:lambdauusize(l,2),l)=[lambdauu(1:lambdauusize(l,1),1:lambdauusize(l,2),l);lambdaru(1:lambdarusize(l,1),1:lambdarusize(l,2),l)];
    U(1:lambdauusize(l,1)+lambdarusize(l,1),lambdauusize(l,2)+psrsize(l+1,2)+1:lambdauusize(l,2)+psrsize(l+1,2)+lambdaussize(l,2),l)=[lambdaus(1:lambdaussize(l,1),1:lambdaussize(l,2),l);lambdars(1:lambdarssize(l,1),1:lambdarssize(l,2),l)];  
    
     L(end-lambdadssize(l,1)+1:end,end-lambdadssize(l,2)-lambdaddsize(l,2)+1:end,l)=[lambdads(1:lambdadssize(l,1),1:lambdadssize(l,2),l),lambdadd(1:lambdaddsize(l,1),1:lambdaddsize(l,2),l)];
end

%there is only an S matrix at the end of omega
S(:,:,end)=[pss(1:psssize(end,1),1:psssize(end,2),end),psd(1:psdsize(end,1),1:psdsize(end,2),end);pszids(1:pszidssize(end-1,1),1:pszidssize(end-1,2),end-1),zeros(pszidssize(end-1,1),psdsize(end,2))];

%creating omega
omega=zeros(statenum*thresnum);
for m=1:thresnum-1
    omega(statenum*(m-1)+1:statenum*m,statenum*(m-1)+1:statenum*m)=S(:,:,m);
    omega(statenum*(m-1)+1:statenum*m,statenum*m+1:statenum*(m+1))=U(:,:,m);
    omega(statenum*m+1:statenum*(m+1),statenum*(m-1)+1:statenum*m)=L(:,:,m);
end


omega(end-statenum+1:end,end-statenum+1:end)=S(:,:,end);

if abs(norm(sum(omega,2)-ones(size(omega,1),1)))>1e-9 
    warning('omega is not stochastic')
end

%checking
%disp('the eigenvalue of omega which is closest to one')
%[sv, se]=eig(omega'); se=diag(se);
%min(abs(se(:)))

%we want to restrict the states according to (42). for this we will order
%the blockmatrices of omega such that it will have the following structure:
%omega=[A, E;
%        F, G] Where A is a block-diagonal matrix which consist of the
%        different pss matrices
I=eye(thresnum*statenum);
meret=0;
for l=1:thresnum
    temp0=I(:,1:meret);
    temp1=I(:,(l-1)*statenum+psziussize(l,1)+pszirssize(l,1)+1:(l-1)*statenum+psziussize(l,1)+pszirssize(l,1)+psssize(l,2));
    temp2=[I(:,meret+1:(l-1)*statenum+psziussize(l,1)+pszirssize(l,1)),I(:,(l-1)*statenum+psziussize(l,1)+pszirssize(l,1)+psssize(l,2)+1:end)];
    I=[temp0,temp1,temp2];
    meret=meret+psssize(l,2);
end

omega=I'*omega*I;


%restricting the state-space
A=omega(1:meret,1:meret);
E=omega(1:meret,meret+1:end);
F=omega(meret+1:end,1:meret);
G=omega(meret+1:end,meret+1:end);

omega=A+E*((eye(size(G))-G)\eye(size(G)))*F;


if abs(norm(sum(omega,2)-ones(size(omega,1),1)))>1e-9 
    warning('the restriced omega is not stochastic')
end

%----------------caclulation of the probabilities--------------------------

% the linear equations we need the diagonal values of Q for the sticky
% states
stickdiag=zeros(size(omega,1),1);
%sticky states at the lowe buffer
index=0;
for l=1:statenum
    if R(l,l,1)<0
        index=index+1;
        stickdiag(index)=-Q(l,l,1);
    end
end
%intermediate
for l=2:thresnum-1
    for m=1:statenum
        if R(m,m,l-1)>0 && R(m,m,l)<0
            index=index+1;
            stickdiag(index)=-Q(m,m,l);
        end
    end
end
%end
for l=1:statenum
    if R(l,l,end)>0
        index=index+1;
        stickdiag(index)=-Q(l,l,end);
    end
end

if index ~= size(omega,1)
    warning('the number of sticky states is not correct for all the pss matrices')
end

FI=diag(stickdiag);

%solvin the linear equations
LINEQ=FI*(omega-eye(size(omega)));

p=null(LINEQ')'; %this vector contains the probabilities at the threshold levels

if size(p,1)==0
    [sv, se]=eig(LINEQ');
    se=abs(diag(se));
    meg=find(se(:)==min(se(:)));
    p=sv(:,meg(1))';
end


accuracy=norm(p*LINEQ);

if norm(abs(p)-p)>1e-9 && norm(abs(p)+p)>1e-9
    disp('there are negative probabilities in the results')
    %this is common when omega has multiple eigenvalues close to 1
    k=0;
    [sv, se]=eig(LINEQ');
    se=abs(diag(se));
    meg=find(se(:)==min(se(:)));
    while norm(abs(p)-p)>1e-9 && norm(abs(p)+p)>1e-9 && k<size(se,1)
    k=k+1;
    se(meg(1))=max(se(:));
    meg=find(se(:)==min(se(:)));
    p=sv(:,meg(1))';
    end
    if k~=size(se,1)
        disp('the accuracy of the next result')
        disp(se(meg(1)))
    else
        disp('the results are not eccurate enough')
    end
end

p=abs(p);

%extracting the data
pmatrix=zeros(thresnum,statenum);
index=0;
for l=1:thresnum
    pmatrix(l,1:psssize(l,1))=p(index+1:index+psssize(l,1));
    index=index+psssize(l,1);
end

%-------calculation of the probability density near thershold levels-------
%we need to create an other huge matrix (equations 32-35)
disp(' ');
disp('calculating the probability density values near threshold levels')



LINEQ=zeros(statenum*thresnum);
VALUE=zeros(1,thresnum*statenum);

%first threshold
cup=[];
cdown=[];
for l=1:statenum
    if R(l,l,1)>0 && R(l,l,2)>0
        cup=[cup,R(l,l,2)];
    elseif R(l,l,1)<0 && R(l,l,2)<0 
        cdown=[cdown,abs(R(l,l,1))];
    end
end
CUP=diag(cup);
CDOWN=diag(cdown);

%taking some values of Q
%at the bottom
TRANS=TRANSmatrix(:,:,1);
upnum=upnumvector(1);
sticknum=sticknumvector(1);
T=TRANS'*Q(:,:,1)*TRANS;
Tsu1=T(upnum+1:end,1:upnum);

%one level higher
TRANS=TRANSmatrix(:,:,2);
upnum=upnumvector(2);
sticknum=sticknumvector(2);
downnum=downnumvector(2);
repulsivenum=repulsivenumvector(2);
T=TRANS'*Q(:,:,2)*TRANS;
Tsu2=T(upnum+1:upnum+sticknum,1:upnum);
Tsd2=T(upnum+1:upnum+sticknum,upnum+sticknum+1:upnum+sticknum+downnum);
Tsr2=T(upnum+1:upnum+sticknum,upnum+sticknum+downnum+1:end);

%first LINEQ
LINEQ(1:size(CUP,1)+size(CDOWN,1),1:size(CUP,2))=[-CDOWN*pszidu(1:pszidusize(1,1),1:pszidusize(1,2),1);CUP];
VALUE(1:size(CUP,1))=pmatrix(1,1:psssize(1,1))*Tsu1*lambdauu(1:lambdauusize(1,1),1:lambdauusize(1,2),1)+pmatrix(2,1:psssize(2,1))*Tsu2;
koord=size(CUP,1);%koord will show how many linear equations were created so far
%second
cdown2=[];
for l=1:statenum
    if R(l,l,2)<0 && R(l,l,3)<0
        cdown2=[cdown2,abs(R(l,l,2))];
    end
end
CDOWN2=diag(cdown2);

LINEQ(1:size(CUP,1)+size(CDOWN,1)+size(CDOWN2,1),koord+1:koord+size(CDOWN,2))=[CDOWN;-CUP*psziud(1:psziudsize(2,1),1:psziudsize(2,2),2);-CDOWN2*lambdadd(1:lambdaddsize(2,1),1:lambdaddsize(2,2),2)];
VALUE(koord+1:koord+size(CDOWN,1))=pmatrix(2,1:psssize(2,1))*(Tsd2+Tsr2*pszird(1:pszirdsize(2,1),1:pszirdsize(2,2),2));
koord=koord+size(CDOWN,1);

%until the thresnum-3nd threshold level
for l=2:thresnum-3
    cup0=[]; cup2=[];
    cdown1=[]; cdown3=[];
    for m=1:statenum
        if R(m,m,l)>0 && R(m,m,l+1)>0
            cup2=[cup2,R(m,m,l+1)];
        elseif R(m,m,l)<0 && R(m,m,l+1)<0
            cdown1=[cdown1,abs(R(m,m,l))];
        end
        if R(m,m,l+1)<0 && R(m,m,l+2)<0
            cdown3=[cdown3,abs(R(m,m,l+1))];
        end
        if R(m,m,l-1)>0 && R(m,m,l)>0
            cup0=[cup0,R(m,m,l)];
        end
    end
    CUP0=diag(cup0); CUP2=diag(cup2);
    CDOWN1=diag(cdown1); CDOWN3=diag(cdown3);
    
    TRANS=TRANSmatrix(:,:,l);
    sticknum=sticknumvector(l);
    repulsivenum=repulsivenumvector(l);
    upnum=upnumvector(l);
    downnum=downnumvector(l);
    Qtemp=TRANS'*Q(:,:,l)*TRANS;
    Tsr1=Qtemp(upnum+1:upnum+sticknum,upnum+sticknum+downnum+1:end);
    
    TRANS=TRANSmatrix(:,:,l+1);
    upnum=upnumvector(l+1);
    sticknum=sticknumvector(l+1);
    downnum=downnumvector(l+1);
    repulsivenum=repulsivenumvector(l+1);
    
    Qtemp=TRANS'*Q(:,:,l+1)*TRANS;
    Tsd2=Qtemp(upnum+1:upnum+sticknum,upnum+sticknum+1:upnum+sticknum+downnum);
    Tsr2=Qtemp(upnum+1:upnum+sticknum,upnum+sticknum+downnum+1:end);
    Tsu2=Qtemp(upnum+1:upnum+sticknum,1:upnum);
    
    %The linear equations:
    %first
    LINEQ(koord+1:koord+size(CDOWN1,1)+size(CDOWN3,1)+size(CUP2,1),koord+1:koord+size(CDOWN1,1))=[CDOWN1;-CUP2*psziud(1:psziudsize(l+1,1),1:psziudsize(l+1,2),l+1);-CDOWN3*lambdadd(1:lambdaddsize(l+1,1),1:lambdaddsize(l+1,2),l+1)];
    VALUE(koord+1:koord+size(CDOWN1,1))=pmatrix(l+1,1:psssize(l+1,1))*Tsd2+pmatrix(l+1,1:psssize(l+1,1))*Tsr2*pszird(1:pszirdsize(l+1,1),1:pszirdsize(l+1,2),l+1);
    
    %second
    LINEQ(koord-size(CUP0,1)+1:koord+size(CDOWN1,1)+size(CUP2,1),koord+size(CDOWN1,1)+1:koord+size(CDOWN1,1)+size(CUP2,1))=[-CUP0*lambdauu(1:lambdauusize(l,1),1:lambdauusize(l,2),l); -CDOWN1*pszidu(1:pszidusize(l,1),1:pszidusize(l,2),l); CUP2];
    VALUE(koord+size(CDOWN1,1)+1:koord+size(CDOWN1,1)+size(CUP2,1))=pmatrix(l+1,1:psssize(l+1,1))*Tsu2+pmatrix(l,1:psssize(l,1))*Tsr1*lambdaru(1:lambdarusize(l,1),1:lambdarusize(l,2),l);
    koord=koord+size(CDOWN1,1)+size(CUP2,1);
end

    %for the end-1 threshold:
    TRANS=TRANSmatrix(:,:,end-1);
    upnum=upnumvector(end-1);
    sticknum=sticknumvector(end-1);
    downnum=downnumvector(end-1);
    repulsivenum=repulsivenumvector(end-1);
       
    Rtemp=TRANS'*abs(R(:,:,end-1))*TRANS;
    CDOWN1=Rtemp(upnum+sticknum+1:upnum+sticknum+downnum,upnum+sticknum+1:upnum+sticknum+downnum);
    
    Rtemp=TRANS'*abs(R(:,:,end))*TRANS;
    CUP2=Rtemp(1:upnum,1:upnum);
    CDOWN2=Rtemp(upnum+sticknum+1:upnum+sticknum+downnum);
    
    Qtemp=TRANS'*Q(:,:,end)*TRANS;
    Tsd1=Qtemp(upnum+1:upnum+sticknum,upnum+sticknum+1:upnum+sticknum+downnum);
    Tsu1=Qtemp(upnum+1:upnum+sticknum,1:upnum);
    
    
    %calculation of Tsd2:
    TRANS=TRANSmatrix(:,:,end);
    upnum=upnumvector(end);
    sticknum=sticknumvector(end);
    downnum=downnumvector(end);
    repulsivenum=repulsivenumvector(end);
    
    Qtemp=TRANS'*Q(:,:,end)*TRANS;
    Tsd2=Qtemp(1:sticknum,sticknum+1:end);
    %first linear equations
    LINEQ(koord+1:koord+size(CDOWN1,1)+size(CUP2,1),koord+1:koord+size(CDOWN1,2))=[CDOWN1;-CUP2*psziud(1:psziudsize(end-1,1),1:psziudsize(end-1,2),end-1)];
    VALUE(koord+1:koord+size(CDOWN1,2))=(pmatrix(end-1,1:psssize(end-1,1))*Tsd1+pmatrix(end,1:psssize(end,1))*Tsd2*lambdadd(1:lambdaddsize(end-1,1),1:lambdaddsize(end-1,2),end-1));
    
    %second linear equation
    TRANS=TRANSmatrix(:,:,end-2);
    upnum=upnumvector(end-2);
    sticknum=sticknumvector(end-2);
    downnum=downnumvector(end-2);
    repulsivenum=repulsivenumvector(end-2);
    
    Qtemp=TRANS'*Q(:,:,end-1)*TRANS;
    Tsr0=Qtemp(upnum+1:upnum+sticknum,upnum+sticknum+downnum+1:end);
    
    Rtemp=TRANS'*abs(R(:,:,end-1))*TRANS;
    CUP0=Rtemp(1:upnum,1:upnum);
 
    LINEQ(koord-size(CUP0,1)+1:koord+size(CUP2,1)+size(CDOWN1,1),koord+size(CDOWN1,2)+1:koord+size(CDOWN1,2)+size(CUP2,2))=[-CUP0*lambdauu(1:lambdauusize(end-2,1),1:lambdauusize(end-2,2),end-2);-CDOWN1*pszidu(1:pszidusize(end-2,1),1:pszidusize(end-2,2),end-2);CUP2];
    VALUE(koord+size(CDOWN1,1)+1:koord+size(CDOWN1,1)+size(CUP2,1))=pmatrix(end-1,1:psssize(end-1,1))*Tsu1+pmatrix(end-2,1:psssize(end-2,1))*Tsr0*lambdaru(1:lambdarusize(end-2,1),1:lambdarusize(end-2,2),end-2);
 
    koord=koord+size(CDOWN1,1)+size(CUP2,1);
    
    LINEQ=LINEQ(1:koord,1:koord);
    VALUE=VALUE(1:koord);
    sol=VALUE/LINEQ;
    
    
    if norm(sol-abs(sol))/norm(sol)>1e-8
       warning('the probability density function has negative values');
    end
    
    %acquiring the probabilities from the results
    pid=zeros(thresnum-2,statenum);
    piu=pid;
    
    koord=0;
    for l=2:thresnum-1
        pid(l-1,1:lambdadssize(l-1,1))=sol(koord+1:koord+lambdadssize(l-1,1));
        koord=koord+lambdadssize(l-1,1);
        piu(l-1,1:lambdauusize(l-1,2))=sol(koord+1:koord+lambdauusize(l-1,2));
        koord=koord+lambdauusize(l-1,2);
        
    end
%---------------calculation of N and everything else-----------------------
disp(' ');
disp('calculation the probability density function')
h=(thres(end)-thres(1))/(n-1);

f=zeros(n,statenum);
koord=0;
level_index=0;
norma=0;%the normalizing constant (numerical normalization)
for x=thres(1):h:thres(end)
    if x>=thres(level_index+1) && x~=thres(end) %if we go over a threshold level we need to recalculate a matrix for N
        level_index=level_index+1;
        cp=[]; cm=[];
        plusvector=zeros(statenum); plusnum=0;
        minususvector=zeros(statenum); minusnum=0;
        for l=1:statenum
            if R(l,l,level_index)>0
                cp=[cp,R(l,l,level_index)];
                plusnum=plusnum+1;
                plusvector(plusnum)=l;
            else
                cm=[cm,-R(l,l,level_index)];
                minusnum=minusnum+1;
                minususvector(minusnum)=l;
            end
        end
        CP=diag(cp); CM=diag(cm);
        
        %we need different parts of Q to calculate K
        TRANS=zeros(statenum);
        for l=1:statenum
            if l<=plusnum
                TRANS(plusvector(l),l)=1;
            else
                TRANS(minususvector(l-plusnum),l)=1;
            end
        end
        RETRANS=TRANS;
        Qtemp=TRANS'*Q(:,:,level_index)*TRANS;
        Tpp=Qtemp(1:plusnum,1:plusnum);
        Tpm=Qtemp(1:plusnum,plusnum+1:end);
        Tmp=Qtemp(plusnum+1:end,1:plusnum);
        Tmm=Qtemp(plusnum+1:end,plusnum+1:end);
    
        %{
         K=CP^(-1)*Tpp+pszipm(1:pszipmsize(level_index,1),1:pszipmsize(level_index,2),level_index)*CM^(-1)*Tmp;
         K2=CM^(-1)*Tmm+pszimp(1:pszimpsize(level_index,1),1:pszimpsize(level_index,2),level_index)*CP^(-1)*Tpm;
         b=thres(level_index+1)-thres(level_index);
    
         MATRIX=eye(statenum);
         MATRIX(1:size(K,1),pszimpsize(level_index,2)+1:end)=expm(K*b)*pszipm(1:pszipmsize(level_index,1),1:pszipmsize(level_index,2),level_index);
         MATRIX(size(K,1)+1:end,1:pszimpsize(level_index,2))=expm(K2*b)*pszimp(1:pszimpsize(level_index,1),1:pszimpsize(level_index,2),level_index);
         MATRIX=MATRIX\eye(statenum);
        %}
        K=CP\Tpp+PSZI(1:PSZIsize(level_index,1),1:PSZIsize(level_index,2),level_index)/CM*Tmp;
         K2=CM\Tmm+PSZI2(1:PSZI2size(level_index,1),1:PSZI2size(level_index,2),level_index)/CP*Tpm;
         b=thres(level_index+1)-thres(level_index);
         
         MATRIX=eye(statenum);
         MATRIX(1:size(K,1),PSZI2size(level_index,2)+1:end)=expm(K*b)*PSZI(1:PSZIsize(level_index,1),1:PSZIsize(level_index,2),level_index);
         MATRIX(size(K,1)+1:end,1:PSZI2size(level_index,2))=expm(K2*b)*PSZI2(1:PSZI2size(level_index,1),1:PSZI2size(level_index,2),level_index);
         MATRIX=MATRIX\eye(statenum);

         %calculatin TRANS to separate the up,repulsive... states
         %for the bottom
         TRANS=TRANSmatrix(:,:,level_index);
         upnum=upnumvector(level_index);
         sticknum=sticknumvector(level_index);
         downnum=downnumvector(level_index);
         repulsivenum=repulsivenumvector(level_index);
                
        TRANS1=TRANS; upnum1=upnum; sticknum1=sticknum; downnum1=downnum; repulsivenum1=repulsivenum;
        
        %Tsu, Tsr and Cu:
        Qtemp=TRANS1'*Q(:,:,level_index)*TRANS1;
        Tsu=Qtemp(upnum1+1:upnum1+sticknum1,1:upnum1);
        Tsr=Qtemp(upnum1+1:upnum1+sticknum1,upnum1+sticknum1+downnum1+1:end);
        Rtemp=TRANS1'*R(:,:,level_index)*TRANS1;
        Cu=Rtemp(1:upnum1,1:upnum1);
        
           
        %up1 and repulsive1 are needed to be separated for Np
        TRANStemp=[TRANS1(:,1:upnum1), TRANS1(:,upnum1+sticknum1+downnum1+1:end)];
        theend=0;
        if size(TRANStemp,2)~=size(TRANStemp,1)
            m=0;
                while theend==0
                  m=m+1;
                  if norm(TRANStemp(m,:))<1e-10
                    TRANStemp=[TRANStemp(1:m-1,:);TRANStemp(m+1:end,:)];
                    m=m-1;
                    if size(TRANStemp,2)==size(TRANStemp,1)
                            theend=1;
                    end                
                  end
                end
        end
        
        TRANS1=TRANStemp;
        
        
        
        %for the top
        TRANS=TRANSmatrix(:,:,level_index+1);
         upnum=upnumvector(level_index+1);
         sticknum=sticknumvector(level_index+1);
         downnum=downnumvector(level_index+1);
         repulsivenum=repulsivenumvector(level_index+1);
        
        TRANS2=TRANS; upnum2=upnum; sticknum2=sticknum; downnum2=downnum; repulsivenum2=repulsivenum;
        
        %Cd, Tsd
        Qtemp=TRANS2'*Q(:,:,level_index)*TRANS2;
        Tsd=Qtemp(upnum2+1:upnum2+sticknum2,upnum2+sticknum2+1:upnum2+sticknum2+downnum2);
        Rtemp=TRANS2'*R(:,:,level_index)*TRANS2;
        Cd=-(Rtemp(upnum2+sticknum2+1:upnum2+sticknum2+downnum2,upnum2+sticknum2+1:upnum2+sticknum2+downnum2));

        
        %repulsive2 and down2 are separated for Nm
        TRANStemp=[TRANS2(:,upnum2+sticknum2+1:upnum2+sticknum2+downnum2), TRANS2(:,upnum2+sticknum2+downnum2+1:end)];
        theend=0;
        if size(TRANStemp,2)~=size(TRANStemp,1)
            m=0;
                while theend==0
                  m=m+1;
                  if norm(TRANStemp(m,:))<1e-10
                    TRANStemp=[TRANStemp(1:m-1,:);TRANStemp(m+1:end,:)];
                    m=m-1;
                    if size(TRANStemp,2)==size(TRANStemp,1)
                            theend=1;
                    end                
                  end
                end
        end
        
        TRANS2=TRANStemp;
        
        
        % normalizing:----------------------------------------------------
        %the analytical normalization is not obvios, we need to separate
        %the eiegenvalue corresponding to 0
        
        a=null(K); b=null(K')'; 
        if isempty(a)
            PI=0;
        else
            PI=a*b/(b*a);
        end
        a=null(K2); b=null(K2')'; 
        if isempty(a)
            PI2=0;
        else
            PI2=a*b/(b*a);
        end
        A=(expm(K*(thres(level_index+1)-thres(level_index)))-eye(size(K))-(thres(level_index+1)-thres(level_index))*PI)/(K-PI);
        B=(expm(K2*(thres(level_index+1)-thres(level_index)))-eye(size(K2))-(thres(level_index+1)-thres(level_index))*PI2)/(K2-PI2);
        
        
        %A=K\expm(K*(thres(level_index+1)-thres(level_index)));
        %B=K2\expm(K2*(thres(level_index+1)-thres(level_index)));
        intN=MATRIX*[A A*PSZI(1:PSZIsize(level_index,1),1:PSZIsize(level_index,2),level_index);
              B*PSZI2(1:PSZI2size(level_index,1),1:PSZI2size(level_index,2),level_index), B]; %the integrand of N
          
        
        Np=intN(1:plusnum,:); Nm=intN(plusnum+1:end,:);
    
        Np=TRANS1'*Np;
        Nu=Np(1:upnum1,:);
        Nr=Np(upnum1+1:end,:);
    
        Nm=TRANS2'*Nm;
        Nd=Nm(1:downnum2,:);
         
    
        %for the normalization we need
        Rtemp=zeros(statenum);
        Rtemp(1:size(CP,1),1:size(CP,1))=CP;
        Rtemp(size(CP,1)+1:size(CP,1)+size(CM,1),size(CP,1)+1:size(CP,1)+size(CM,1))=CM;
          
        if level_index ==1
            norma=norma+sum((pmatrix(1,1:psssize(1,1))*Tsu*Nu+pid(1,1:lambdadssize(1,1))*Cd*Nd)/Rtemp); %the RETRANS transforms back the states to the original order        
        elseif level_index<=thresnum-2
            norma=norma+sum((piu(level_index-1,1:lambdauusize(level_index-1,2))*Cu*Nu+pmatrix(level_index,1:psssize(level_index,1))*Tsr*Nr+pid(level_index,1:lambdadssize(level_index,1))*Cd*Nd)/Rtemp);
        else% in the last domain
            norma=norma+sum((piu(level_index-1,1:lambdauusize(level_index-1,2))*Cu*Nu+pmatrix(level_index,1:psssize(level_index,1))*Tsr*Nr+pmatrix(end,1:psssize(end,1))*Tsd*Nd)/Rtemp);
        end
        %-----------------------------------------------------------------
    end
    %Calculation of N(x) for different values of x
    %calculates A and B to slightly speed up the calculation of N
    A=expm(K*(x-thres(level_index)));
    B=expm(K2*(thres(level_index+1)-x));

    N=MATRIX*[A A*PSZI(1:PSZIsize(level_index,1),1:PSZIsize(level_index,2),level_index);
              B*PSZI2(1:PSZI2size(level_index,1),1:PSZI2size(level_index,2),level_index), B];
          
    Np=N(1:plusnum,:); Nm=N(plusnum+1:end,:);
    Np=TRANS1'*Np;
    Nu=Np(1:upnum1,:);
    Nr=Np(upnum1+1:end,:);
    
    Nm=TRANS2'*Nm;
    Nd=Nm(1:downnum2,:);

    %the values of the function
    koord=koord+1;
    if x==thres(level_index)
        f(koord,:)=NaN;
    else
        if level_index ==1
        f(koord,:)=(pmatrix(1,1:psssize(1,1))*Tsu*Nu+pid(1,1:lambdadssize(1,1))*Cd*Nd)/Rtemp*RETRANS'; %the RETRANS transforms back the states to the original order        
        elseif level_index<=thresnum-2
       f(koord,:)=(piu(level_index-1,1:lambdauusize(level_index-1,2))*Cu*Nu+pmatrix(level_index,1:psssize(level_index,1))*Tsr*Nr+pid(level_index,1:lambdadssize(level_index,1))*Cd*Nd)/Rtemp*RETRANS';
        else% in the last domain
        f(koord,:)=(piu(level_index-1,1:lambdauusize(level_index-1,2))*Cu*Nu+pmatrix(level_index,1:psssize(level_index,1))*Tsr*Nr+pmatrix(end,1:psssize(end,1))*Tsd*Nd)/Rtemp*RETRANS';
        
        end
    end    
   
end

%numerical normalization
norma=norma+sum(sum(pmatrix));
f=f/norma;
pmatrix=pmatrix/norma;


            
%-----------------------post processing------------------------------------
x=thres(1):h:thres(end);


disp(' ');
disp('Accuracy:')
disp(accuracy)
disp('Probabilities at threshold levels:')
disp('threshold     probability')

if thresnumtemp ~= 0 %if the original thresnum was ==2 or ==3
    thres=thres(2:end-1);
    pmatrix=[pmatrix(1,:);pmatrix(3:end-2,:);pmatrix(end,:)];
else
    thresnumtemp=thresnum;
end
%reordering pmatrix back to the original order of the states
pmatrix2=zeros(size(pmatrix));
for l=1:thresnumtemp
    index=1;
    for m=1:statenum
        if l==1
            if Rorig(m,m,l)<0
                pmatrix2(l,m)=pmatrix(l,index);
                index=index+1;
            end
        elseif l<thresnumtemp
            if Rorig(m,m,l)<0 && Rorig(m,m,l-1)>0
                pmatrix2(l,m)=pmatrix(l,index);
                index=index+1;
            end
        else
            if Rorig(m,m,l-1)>0
                pmatrix2(l,m)=pmatrix(l,index);
                index=index+1;
            end
        end
    end
end
pmatrix=pmatrix2;

disp([thres(:),sum(pmatrix')'])


