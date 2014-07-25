%  [B, n] = MStaircase(X, Z, precision)
%  
%  Computes a smaller representation using the staircase 
%  algorithm.
%  
%  Notes
%  -----
%  This function should not be called directly.
%  It is used by 'MinimalRepFromME' and 'MinimalRepFromRAP'.
%  
%  References
%  ----------
%  .. [1]  P. Buchholz, M. Telek, "On minimal representation 
%          of rational arrival processes." Madrid Conference
%          on Qeueuing theory (MCQT), June 2010.

function [B,n] = MStaircase (X,Z, precision)

    if nargin<3
        precision = 1e-8;
    end

    msize = length(X);
    m = size(X{1},1);
    U = eye(m);

    
    ranksum=0; %The sum of the ranks calculated in every loop
    crit=1;    %The stopping criteria
    while crit
        r=rank(Z,precision);
        ranksum=ranksum+r;
        [Ui,S,T] = svd(Z);

        %Calculation of the new U,X,Z matrices
        Transf=eye(ranksum-r+size(Ui,1));
        Transf(end-size(Ui,1)+1:end,end-size(Ui,1)+1:end)=Ui';
        U=(Transf*U')';

        for ii=1:msize
            TEMP = Ui'*X{ii}*Ui;
            X{ii} = TEMP(r+1:end,r+1:end);
            if ii==1
                Z = TEMP(r+1:end,1:r);
            else
                Z = horzcat(Z,TEMP(r+1:end,1:r));
            end
        end

       if norm(Z)<precision || rank(Z,precision)==m-ranksum
            crit=0;
        end
    end
    
    n=ranksum;
    I=eye(m,m);
    if norm(Z)<precision 
        n=ranksum;
        x=sum(U',2);
        x=x(1:n);

        %does x have a 0 value somewhere
        yes=0;
        zeroloc=zeros(n,1);
        nonzero=0; %this will indicate a row of x for which x's value is non-zero
        for l=1:n
           if abs(x(l))<precision
               yes=1;
               zeroloc(l,1)=1;
           elseif nonzero==0
               nonzero=l;
           end
        end

        R=eye(n,n);
        if yes
            for l=1:n;
               if zeroloc(l,1)==1;
                R(l,nonzero)=1;
               end
            end
        end
        
        y=R*x;
        Gamma=diag(y);
        TEMP1=I;
        TEMP1(1:n,1:n)=inv(Gamma);
        TEMP2=I; 
        TEMP2(1:n,1:n)=R;
        B=inv(TEMP1*TEMP2*U');
    elseif rank(Z,precision)==m-ranksum
        B=I;
        n=m;
    end
end
