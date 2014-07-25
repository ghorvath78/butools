function [res] = ordereigs(A)
%The output of ORDEREIGS is the eigenvalues of a matrix in Schur-form

matrixsize=size(A,1);
res=zeros(matrixsize,1);

l=1;
while l<=matrixsize
    if l==matrixsize
        res(l)=A(l,l);
        l=l+1;
    elseif abs(A(l+1,l))<1e-9*(abs(A(l,l))+abs(A(l+1,l+1)))
        res(l)=A(l,l);
        l=l+1;
    else
        TEMP=A(l:l+1,l:l+1);
        res(l:l+1)=eig(TEMP);
        l=l+2;
    end
end


end

