%  [G0, G1] = CanonicalFromDMAP2(D0, D1, prec)
%  
%  Returns the canonical form of an order-2 discrete Markovian
%  arrival process.
%  
%  Parameters
%  ----------
%  D0 : matrix, shape (2,2)
%      The D0 matrix of the DMAP(2)
%  D1 : matrix, shape (2,2)
%      The D1 matrix of the DMAP(2)
%  prec : double, optional
%      Numerical precision to check the input, default 
%      value is 1e-14
%  
%  Returns
%  -------
%  G0 : matrix, shape (1,2)
%      The D0 matrix of the canonical DMAP(2)
%  G1 : matrix, shape (2,2)
%      The D1 matrix of the canonical DMAP(2)

function [g0,g1]=CanonicalFromDMAP2(d0, d1)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckDMAPRepresentation(d0,d1)
        error('DMAP2Canonical: Input isn''t a valid DMAP representation!');
    end

    if size(d0,1)~=2
        error('DMAP2Canonical: Size is not 2!');
    end

    ev=EigSort(eig(d0));
    s1=ev(1);
    s2=ev(2);

    if s2 >= 0
        [g0,g1]=CanonicalFromMAP2(d0-eye(2),d1);
        g0=g0+eye(2);
        return;
    end

    %s2 is negative

    av=DRPSolve(inv(eye(2)-d0)*d1);
    gamma=EigSort(eig(inv(eye(2)-d0)*d1));
    gamma=gamma(2);

    w1=1/(s1-s2)*(sum(d0,2)-s2*[1;1]);
    w2=[1;1]-w1;

    W=[w1 w2];
    A=(1-s1)*(av*W);
    a1=A(1);

    if gamma>=0
        a=-(1/(2*(-1+s1)*(-1+s1+s2)^2))*(1-4*s1+a1*s1+5*s1^2-a1*s1^2-2*s1^3-2*s2-a1*s2+5*s1*s2-3*s1^2*s2+s2^2+a1*s2^2-s1*s2^2-gamma+3*s1*gamma-a1*s1*gamma-3*s1^2*gamma+a1*s1^2*gamma+s1^3*gamma+s2*gamma+a1*s2*gamma-2*s1*s2*gamma+s1^2*s2*gamma-a1*s2^2*gamma+sqrt((-1+s1+s2)^2*((-1+s1^2*(-2+gamma)+gamma+s2*(1+a1-a1*gamma)+s1*(3-a1-s2-2*gamma+a1*gamma))^2-4*(-1+s1)*(-s1^3*(-1+gamma)+a1*(-1+s2)*s2*(-1+gamma)+s1^2*(-2+a1+s2+2*gamma-a1*gamma)+s1*(1-a1-s2-gamma+a1*gamma)))));
        b=1+(a*(-1+s1+s2-s1*s2)*gamma)/((a-1)*(-s1*s2+a*(-1+s1+s2)));

        g0=[s1+s2 a*(1-s1-s2); s1*s2/(a*(s1+s2-1)) 0];
        g1=[(1-a)*(1-s1-s2) 0; b*(1+s1*s2/(a*(1-s1-s2))) (1-b)*(1+s1*s2/(a*(1-s1-s2)))];
    else        
        %gamma<0
        a=(a1*s1-a1*s1^2+s2-a1*s2-3*s1*s2+2*s1^2*s2-s2^2+a1*s2^2+s1*s2^2+s1*gamma-a1*s1*gamma-2*s1^2*gamma+a1*s1^2*gamma+s1^3*gamma+a1*s2*gamma-a1*s2^2*gamma+sqrt(-4*(-1+s1)*s1*s2*(-1+s1+s2)*(a1*(s1-s2)*(-1+gamma)+(-1+s1)*(s2+(-1+s1)*gamma))+(a1*(-s1+s1^2+s2-s2^2)*(-1+gamma)+(-1+s1)*((-1+2*s1)*s2+s2^2+(-1+s1)*s1*gamma))^2))/(2*(-1+s1+s2)*(a1*(s1-s2)*(-1+gamma)+(-1+s1)*(s2+(-1+s1)*gamma)));
        b=-((a*(1-s1)*(1-s2)*gamma)/((a-1)*(-a+a*s1+a*s2-s1*s2)));

        g0=[s1+s2 a*(1-s1-s2); s1*s2/(a*(s1+s2-1)) 0];
        g1=[0 (1-a)*(1-s1-s2); b*(1-s1*s2/(a*(s1+s2-1))) (1-b)*(1-s1*s2/(a*(s1+s2-1)))];
    end
end