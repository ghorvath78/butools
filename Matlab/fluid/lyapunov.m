function X = lyapunov(A, B, C)
%LYAP solves the equation A*X + X*B + C = 0 for X based on the
%algorithm in Danny C. Sorensen and Yunkai Zhou - Direct methods for matrix
%Sylvester and Lyapunov equations
[Q1,R1]=schur(A);
[Q2,R2]=schur(B);

D=Q1'*C*Q2;
Rsq=R1*R1;
I=eye(size(A));

Xtemp=zeros(size(R1,2),size(R2,1));

l=1;
while l<size(B,1)+1
    if  l==size(B,1) || R2(l+1,l)<1e-9*max(abs([R2(l,l),R2(l+1,l+1)]))
        b=-D(:,l)-Xtemp*R2(:,l);
        LINEQ=R1+R2(l,l)*I;
        Xtemp(:,l)=LINEQ\b;
        l=l+1;
    else
        r11=R2(l,l); r12=R2(l,l+1);
        r21=R2(l+1,l); r22=R2(l+1,l+1);
        b=-D(:,l:l+1)-Xtemp(:,1:l-1)*R2(1:l-1,l:l+1);
        b=[R1*b(:,1)+r22*b(:,1)-r21*b(:,2), R1*b(:,2)+r11*b(:,2)-r12*b(:,1)];
        LINEQ=Rsq+(r11+r22)*R1+(r11*r22-r12*r21)*I;
        Xtemp(:,l:l+1)=LINEQ\b;
        l=l+2;
    end
end

X=Q1*Xtemp*Q2';

end

