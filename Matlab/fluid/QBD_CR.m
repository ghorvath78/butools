function [G,R,U]=QBD_CR(A0,A1,A2)
%QBD_CR Solves the equation G=A0+A1*G+A2*G^2, where the row sums of the
%matrix A0+A1+A2 is zero

A1k=A1-eye(size(A1));
A1v=A1k;
A2v=A2;
A0v=A0;

stop_criteria=0;
while stop_criteria == 0
    A1ktemp=A1k-A2v/A1v*A0v;
    A1vtemp=A1v-A2v/A1v*A0v-A0v/A1v*A2v;
    A2vtemp=-A2v/A1v*A2v;
    A0vtemp=-A0v/A1v*A0v;

    if norm(A1ktemp-A1k)/norm(A1k)<1e-11
        stop_criteria=1;
    end
    A1k=A1ktemp; A1v=A1vtemp; A2v=A2vtemp; A0v=A0vtemp;
end

G=-A1k\A0;
R=A2*(eye(size(G))-A1-A2*G)^(-1);
U=A1+A2*G;

end

