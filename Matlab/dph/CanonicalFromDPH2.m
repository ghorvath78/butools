%  [beta, B] = CanonicalFromDPH2(alpha, A, prec)
%  
%  Returns the canonical form of an order-2 discrete phase-type 
%  distribution.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,2)
%      Initial vector of the discrete phase-type distribution
%  A : matrix, shape (2,2)
%      Transition probability matrix of the discrete phase-type
%      distribution
%  prec : double, optional
%    Numerical precision for checking the input, default value
%    is 1e-14
%  
%  Returns
%  -------
%  beta : matrix, shape (1,2)
%    The initial probability vector of the canonical form
%  B : matrix, shape (2,2)
%    Transition probability matrix of the canonical form

function [beta,B] = CanonicalFromDPH2 (alpha,A)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMGRepresentation(alpha,A)
        error('CanonicalFromDPH2: Input isn''t a valid MG distribution!');
    end

    if BuToolsCheckInput && (size(A,1)~=2 || size(A,2)~=2)
        error('CanonicalFromDPH2: Dimension must be 2!');
    end

    lambda=EigSort(eig(A));
    e=[1;1];
    p1=alpha*(e-A*e);

    if lambda(1)>0 && lambda(2)>0 && lambda(1) ~= lambda(2)
        d1=(1-lambda(1))*(1-p1-lambda(2))/(lambda(1)-lambda(2));
        d2=p1-d1;
	    beta=[d1*(lambda(1)-lambda(2))/((1-lambda(1))*(1-lambda(2))),(d1+d2)/(1-lambda(2))];
        B=[lambda(1),1-lambda(1);0,lambda(2)];
    elseif lambda(1)>0 && lambda(1)==lambda(2)
        d2=p1;
        d1=(1-lambda(1))*(1-d2-lambda(1))/lambda(1);
	    beta=[d1*lambda(1)/(1-lambda(1))^2,d2/(1-lambda(1))];
        B=[lambda(1),1-lambda(1);0,lambda(1)];
    elseif lambda(1)>0
        d1=(1-lambda(1))*(1-p1-lambda(2))/(lambda(1)-lambda(2));
        d2=p1-d1;
	    beta=[(d1*lambda(1)+d2*lambda(2))/((1-lambda(1))*(1-lambda(2))),(d1+d2)*(1-lambda(1)-lambda(2))/((1-lambda(1))*(1-lambda(2)))];
        B=[lambda(1)+lambda(2),1-lambda(1)-lambda(2);lambda(1)*lambda(2)/(lambda(1)+lambda(2)-1),0];
    end
end

